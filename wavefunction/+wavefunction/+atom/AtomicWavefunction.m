classdef AtomicWavefunction < handle
    %ATOMICWAVEFUNCTION Summary of this class goes here

    properties (Access = public)
        lMax (1, 1) uint16 = 0
        gk (:, 1) Ggrid
    end

    properties (Access = public)
        crystal Crystal
        Ylm (:, :) double
        PpData struct
        chiq (:, 1) cell
    end

    methods
        function this = AtomicWavefunction(crystal)
            %ATOMICWAVEFUNCTION Construct an instance of this class
            arguments
                crystal Crystal {mustBeNonempty}
            end

            this.crystal = crystal;
            this.computelMax(crystal);
            % 计算 k+G 得到的网格
            this.gk = Ggrid(crystal).addKPoints(crystal.kpts);

            this.computeYlm();
            this.interpolateChi();
        end
    end

    methods (Access = private)
        function computelMax(this, crystal)
            arguments
                this
                crystal Crystal {mustBeNonempty}
            end

            this.PpData = struct();

            symbol = {crystal.atoms.symbol};
            for i = 1:length(symbol)
                pp = pseudopotential.PpData(symbol{i});
                if pp.info.number_of_wfc == 0
                    error('ProcessSuite:WaveFunction:MissingPSWFC', ...
                        'UPF file for element %s must contain PP_PSWFC data.', symbol{i});
                end
                this.lMax = max(this.lMax, max(pp.pswfc.l));
                this.PpData.(symbol{i}) = pp;
            end
        end

        function Ylm = computeYlm(this)
            arguments
                this
            end

            if this.lMax == 0
                Ylm = ones(this.gk.ng, 1) * sqrt(1/(4*pi));
                return
            else
                Ylm = zeros(this.gk.ng, (this.lMax+1)^2, 'double');
            end

            % Compute cos(θ) and sin(θ)
            sqrtGkk = sqrt(this.gk.gkk);
            cosTheta = zeros(size(sqrtGkk));
            mask = sqrtGkk >= eps;
            cosTheta(mask) = this.gk.gkz(mask) ./ sqrtGkk(mask);
            sinTheta = sqrt(max(0, 1-cosTheta.^2));

            % Compute φ for later use in cos(mφ) and sin(mφ)
            condition1 = this.gk.gkx > eps;
            condition2 = this.gk.gkx < -eps;
            condition3 = ~condition1 & ~condition2;

            phi = zeros(this.gk.ng, 1);
            phi(condition1) = atan(this.gk.gky(condition1) ./ this.gk.gkx(condition1));
            phi(condition2) = atan(this.gk.gky(condition2) ./ this.gk.gkx(condition2)) + pi;
            phi(condition3) = pi/2 .* sign(this.gk.gky(condition3));

            % Init the first few values of Qlm
            Qlm = zeros(this.gk.ng, (this.lMax+1)^2, 'double');
            Qlm(:, 1) = sqrt(1.0/4*pi); % l = 0, m = 0
            Qlm(:, 2) = sqrt(3.0/4*pi) * cosTheta; % l = 1, m = 0
            Qlm(:, 3) = -sqrt(3.0/8*pi) * sinTheta; % l = 1, m = 1

            lmax = this.lMax;
            for ig = 1:this.gk.ng
                % Init the first few values of Ylm
                Ylm(ig, 1) = Qlm(ig, 1); % l = 0, m = 0
                Ylm(ig, 2) = Qlm(ig, 2); % l = 1, m = 0
                Ylm(ig, 3) = sqrt(2.0) * Qlm(ig, 3) * cos(phi(ig)); % l = 1, m = 1

                Qlm(ig, 4) = sqrt(3.0/8*pi) * sinTheta(ig); % l = 1, m = -1
                Ylm(ig, 4) = sqrt(2.0) * Qlm(ig, 3) * sin(phi(ig)); % l = 1, m = -1

                for l = 2:lmax
                    for m = 0:l-2
                        lm  = l^2 + 2*m; % (l, m)
                        lm1 = (l-1)^2 + 2*m; % (l-1, m)
                        lm2 = (l-2)^2 + 2*m; % (l-2, m)

                        % Compute the Qlm value using the recursion formula for 0 <= m <= l-2
                        Qlm(ig, lm) = sqrt((4*l^2-1)/(l^2 - m^2)) * cosTheta(ig) * Qlm(ig, lm1) ...
                            - sqrt(((2*l+1)*(l+m-1)*(l-m-1))/((2*l-3)*(l^2-m^2))) * Qlm(ig, lm2);

                        if m == 0
                            % Update Ylm for (l, 0)
                            Ylm(ig, lm) = Qlm(ig, lm);
                        else
                            % Update Ylm for (l, m) and (l, -m)
                            Ylm(ig, lm) = sqrt(2.0) * Qlm(ig, lm) * cos(m * phi(ig)); % Y(l, m)
                            Ylm(ig, lm+1) = sqrt(2.0) * Qlm(ig, lm) * sin(m * phi(ig)); % Y(l, -m)
                        end
                    end

                    % Handle the special case when m = l and m = l-1
                    lm  = l^2 + 2*l; % (l, l)
                    lm1 = l^2 + 2*(l-1); % (l, l-1)
                    lm2 = (l-1)^2 + 2*(l-1); % (l-1, l-1)

                    % Update Qlm and Ylm for m = l
                    Qlm(ig, lm) = -sqrt((2*l+1)/(2*l)) * sinTheta(ig) * Qlm(ig, lm2);
                    Ylm(ig, lm) = sqrt(2.0) * Qlm(ig, lm) * cos(m * phi(ig)); % Y(l, l)
                    Ylm(ig, lm+1) = sqrt(2.0) * Qlm(ig, lm) * cos(m * phi(ig)); % Y(l, -l)
                    % Update Qlm and Ylm for m = l-1
                    Qlm(ig, lm1)  = sqrt(2*l+1) * cosTheta(ig) * Qlm(ig, lm2);
                    Ylm(ig, lm1) = sqrt(2.0) * Qlm(ig, lm1) * cos(m * phi(ig)); % Y(l, l-1)
                    Ylm(ig, lm1+1) = sqrt(2.0) * Qlm(ig, lm1) * cos(m * phi(ig)); % Y(l, -(l-1))
                end
            end
            this.Ylm = Ylm;
        end

        function interpolateChi(this)
            arguments
                this
            end

            % Define spherical Bessel function j_v(z)
            function jv = sphericalBesselJ(v, z)
                jv = zeros(size(z));

                small = abs(z) < eps;
                notSmall = ~small;

                % Handle z >= eps using the standard formula
                jv(notSmall) = sqrt(pi/2) ./ sqrt(z(notSmall)) .* besselj(v+0.5, z(notSmall));

                % Handle z < eps based on the value of v
                if any(small)
                    switch v
                        case -1
                            error('sphericalBesselJ: Undefined for v = -1 when z < eps.');
                        case 0
                            jv(small) = 1;
                        otherwise
                            jv(small) = 0;
                    end
                end
            end

            atomTypes = fieldnames(this.PpData);
            this.chiq = cell(length(atomTypes), 1);
            for i = 1:length(atomTypes)
                atomType = atomTypes{i};

                % Generate q grid (0:dq:max(|q+G|)) with dq=0.01
                dq = 0.01;
                qgSqrt = sqrt(this.gk.gkk);
                qGrid = 0:dq:(max(qgSqrt)+4*dq);

                % Simpson integrate and Lagrange interpolation
                for index = 1:this.PpData.(atomType).info.number_of_wfc
                    l = this.PpData.(atomType).pswfc.l(index);

                    % Limit the radial grid in the psedopotentials files up to 10 a.u.
                    % to cut off the numerical noise arising from the large-r tail
                    % in cases like the integration of V_loc-Z/r
                    r = this.PpData.(atomType).r(this.PpData.(atomType).r <= 10.0);
                    rab = this.PpData.(atomType).rab(1:length(r));

                    % Simpson integrate to compute the table for interpolation
                    result = zeros(length(qGrid), length(r));
                    for q = 1:length(qGrid)
                        fun = sphericalBesselJ(l, qGrid(q) * r) .* ...
                            this.PpData.(atomType).pswfc.chi(1:length(r), index) .* r;
                        if mod(length(r), 2) == 1
                            factor = [1, abs(mod(2:(length(r)-1), 2) - 2) * 2, 1];
                            range = 1:length(r);
                        else
                            factor = [1, abs(mod(2:(length(r)-1), 2) - 2) * 2, -1];
                            range = [1:(length(r)-1), length(r)-1];
                        end
                        result(q, :) = 1.0/3.0 * sum(factor .* fun(range) .* rab(range));
                    end

                    % Final interpolation table
                    resultTable = result * (4*pi/sqrt(this.crystal.vol));

                    % Lagrange interpolation for all |q+G|
                    px = qgSqrt ./ dq - floor(qgSqrt ./ dq);
                    ux = 1.0 - px;
                    vx = 2.0 - px;
                    wx = 3.0 - px;

                    i0 = floor(qgSqrt ./ dq) + 1;
                    i1 = i0 + 1;
                    i2 = i0 + 2;
                    i3 = i0 + 3;

                    this.chiq{i, 1} = ...
                        resultTable(i0, :) .* ux .* vx .* wx ./ 6.0 + ...
                        resultTable(i1, :) .* px .* vx .* wx ./ 2.0 - ...
                        resultTable(i2, :) .* px .* ux .* wx ./ 2.0 + ...
                        resultTable(i3, :) .* px .* ux .* vx ./ 6.0;
                end
            end
        end
    end
end

