classdef AtomicWavefunction < handle
    %ATOMICWAVEFUNCTION Summary of this class goes here

    properties (Access = private)
        lMax (1, 1) uint16 = 0
        gk (:, 1) Ggrid
        Ylm (:, :) double
    end

    methods
        function this = AtomicWavefunction(crystal)
            %ATOMICWAVEFUNCTION Construct an instance of this class
            arguments
                crystal Crystal {mustBeNonempty}
            end

            this.computelMax(crystal);
            % 计算 k+G 得到的网格
            this.gk = Ggrid(crystal).addKPoints(crystal.kpts);

            this.Ylm = this.computeYlm();
        end
    end

    methods (Access = private)
        function computelMax(this, crystal)
            arguments
                this
                crystal Crystal {mustBeNonempty}
            end

            symbol = {crystal.atoms.symbol};
            for i = 1:length(symbol)
                pp = pseudopotential.PpData(symbol{i});
                if pp.info.number_of_wfc == 0
                    error('ProcessSuite:WaveFunction:MissingPSWFC', ...
                        'UPF file for element %s must contain PP_PSWFC data.', symbol{i});
                end
                this.lMax = max(this.lMax, max(pp.pswfc.l));
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
            if sqrtGkk < eps
                cosTheta = 0;
            else
                cosTheta = this.gk.gkz ./ sqrtGkk;
            end
            sinTheta = sqrt(max(0, 1-cosTheta.^2));

            % Compute φ for later use in cos(mφ) and sin(mφ)
            condition1 = this.gk.gkx > eps;
            condition2 = this.gk.gkx < -eps;
            condition3 = ~condition1 & ~condition2;

            phi = zeros(this.gk.ng, 1);
            phi(condition1) = atan(this.gk.gky(condition1) ./ this.gk.gkx(condition1));
            phi(condition2) = atan(this.gk.gky(condition2) ./ this.gk.gkx(condition2)) + pi;
            phi(condition3) = sign(pi / 2, this.gk.gky(condition3));

            % Init the first few values of Qlm
            Qlm = zeros(this.gk.ng, (this.lMax+1)^2, 'double');
            Qlm(:, 1) = sqrt(1.0/4*pi); % l = 0, m = 0
            Qlm(:, 2) = sqrt(3.0/4*pi) * cosTheta; % l = 1, m = 0
            Qlm(:, 4) = -sqrt(3.0/8*pi) * sinTheta; % l = 1, m = 1

            lmax = this.lMax;
            for ig = 1:this.gk.ng
                % Init the first few values of Ylm
                Ylm(ig, 1) = Qlm(ig, 1); % l = 0, m = 0
                Ylm(ig, 2) = Qlm(ig, 2); % l = 1, m = 0
                Ylm(ig, 4) = sqrt(2.0) * Qlm(ig, 4) * cos(phi(ig)); % l = 1, m = 1

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
        end
    end
end

