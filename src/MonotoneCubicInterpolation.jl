module MonotoneCubicInterpolation
    export interpolateCubicHermite, interpolate_derivative_CubicHermite
    """
    Interpolate using cubic Hermite splines. The breakpoints in arrays xbp and ybp are assumed to be sorted.
    Evaluate the function in all points of the array xeval.
    Methods:
        "Linear"                yuck
        "FiniteDifference"      classic cubic interpolation, no tension parameter
        "Cardinal"              cubic cardinal splines, uses tension parameter which must be between [0,1]
        "FritschCarlson"        monotonic - tangents are first initialized, then adjusted if they are not monotonic
        "FritschButland"        monotonic - faster algorithm (only requires one pass) but somewhat higher apparent "tension"
        "Steffen"               monotonic - also only one pass, results usually between FritschCarlson and FritschButland
    Sources:
        Fritsch & Carlson (1980), "Monotone Piecewise Cubic Interpolation", doi:10.1137/0717021.
        Fritsch & Butland (1984), "A Method for Constructing Local Monotone Piecewise Cubic Interpolants", doi:10.1137/0905021.
        Steffen (1990), "A Simple Method for Monotonic Interpolation in One Dimension", http:#adsabs.harvard.edu/abs/1990A%26A...239..443S
    """
    function interpolateCubicHermite(xeval, xbp, ybp, method, tension)
        method = lowercase(method)
        # first we need to determine tangents (m)
        n = length(xbp)
        m, delta = calcTangents(xbp, ybp, method, tension)

        c = zeros(n-1)
        d = zeros(n-1)
        for k = 1:n-1
            if method == "linear"
                m[k] = delta[k]
                c[k] = 0
                d[k] = 0
                continue
            end
            xdiff = xbp[k+1] - xbp[k]
            c[k] = (3*delta[k] - 2*m[k] - m[k+1]) / xdiff
            d[k] = (m[k] + m[k+1] - 2*delta[k]) / xdiff / xdiff
        end
        len = length(xeval)
        f = zeros(len)
        k = 1
        for i = 1:len
            x = xeval[i]
            if (x < xbp[1]) || (x > xbp[n])
                throw(string("interpolateCubicHermite: x value ", x ," outside breakpoint range [", xbp[1] ,", " ,xbp[n], "]"))
            end
            while (k < n) && (x > xbp[k+1])
                k+=1
            end
            xdiff = x - xbp[k]
            f[i] = ybp[k] + m[k]*xdiff + c[k]*xdiff*xdiff + d[k]*xdiff*xdiff*xdiff
        end
        return f
    end

    function calcTangents(x, y, method, tension)
        n = length(x)
        delta = zeros(n-1)
        m = zeros(n)
        for k = 1:n-1 #   (var k=0; k < n-1; k++) {
            deltak = (y[k+1] - y[k]) / (x[k+1] - x[k])
            delta[k] = deltak
            if k == 1    # left endpoint, same for all methods
                m[k] = deltak
            else
                if method == "cardinal"
                    m[k] = (1 - tension) * (y[k+1] - y[k-1]) / (x[k+1] - x[k-1])
                elseif method == "fritschbutland"
                    alpha = (1 + (x[k+1] - x[k]) / (x[k+1] - x[k-1])) / 3  # Not the same alpha as below.
                    if delta[k-1] * deltak <= 0
                        m[k] = 0
                    else
                        m[k] = delta[k-1] * deltak / (alpha*deltak + (1-alpha)*delta[k-1])
                    end
                elseif method == "fritschcarlson"
                    # If any consecutive secant lines change sign (i.e. curve changes direction), initialize the tangent to zero.
                    # This is needed to make the interpolation monotonic. Otherwise set tangent to the average of the secants.
                    if delta[k-1] * deltak < 0
                        m[k] = 0
                    else
                        m[k] = (delta[k-1] + deltak) / 2
                    end
                #elseif method == "steffen"
                #    p = ((x[k+1] - x[k]) * delta[k-1] + (x[k] - x[k-1]) * deltak) / (x[k+1] - x[k-1])
                #    m[k] = (sign(delta[k-1]) + sign(deltak)) * min(abs(delta[k-1]), abs(deltak), 0.5*abs(p))
                else    # FiniteDifference
                    m[k] = (delta[k-1] + deltak) / 2
                end
            end
        end
        m[n] = delta[n-1]

        if method != "fritschcarlson"
            return m, delta
        end
        """
        Fritsch & Carlson derived necessary and sufficient conditions for monotonicity in their 1980 paper. Splines will be
        monotonic if all tangents are in a certain region of the alpha-beta plane, with alpha and beta as defined below.
        A robust choice is to put alpha & beta within a circle around origo with radius 3. The FritschCarlson algorithm
        makes simple initial estimates of tangents and then does another pass over data points to move any outlier tangents
        into the monotonic region. FritschButland & Steffen algorithms make more elaborate first estimates of tangents that
        are guaranteed to lie in the monotonic region, so no second pass is necessary.
        """
        # Second pass of FritschCarlson: adjust any non-monotonic tangents.
        for k = 1:n-1 #(var k=0; k < n-1; k++) {
            deltak = delta[k]
            if (deltak == 0)
                m[k] = 0
                m[k+1] = 0
                continue
            end
            alpha = m[k] / deltak
            beta = m[k+1] / deltak
            tau = 3 / sqrt(alpha.^2 + beta.^2);
            if (tau < 1)       # if we're outside the circle with radius 3 then move onto the circle
                m[k] = tau * alpha * deltak
                m[k+1] = tau * beta * deltak
            end
        end
        return m, delta
    end

end # module
