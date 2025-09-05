classdef Projector < handle
    %PROJECTOR Summary of this class goes here

    properties
        originalAtomicWavefunction
    end

    methods
        function this = Projector(crystal)
            %PROJECTOR Construct an instance of this class
            arguments
                crystal Crystal {mustBeNonempty}
            end

            this.originalAtomicWavefunction = wavefunction.atom.AtomicWavefunction(crystal).atomicWavefunction;
        end

        function c = project(this, X)
            arguments
                this 
                X {mustBeA(X, {'Wavefun', 'BlochWavefun'})}
            end

            % Calculate the overlap matrix O
            overlapMatrix = this.originalAtomicWavefunction' * this.originalAtomicWavefunction;

            % Calculate eigenvalues and eigenvectors of the overlap matrix O
            % eigenVectors: columns are eigenvectors of O
            % eigenValues: diagonal matrix containing eigenvalues of O
            [eigenVectors, eigenValues] = eig(overlapMatrix);

            % Calculate inverse square root of eigenvalue matrix
            % This gives O^(-1/2) in diagonal form
            eigenValuesSqrtInv = sqrt(inv(eigenValues));

            % Construct the symmetric orthogonalization matrix O^(-1/2)
            % O^(-1/2) = L * λ^(-1/2) * L^†, where L contains eigenvectors as columns
            symmetricOrthogMatrix = eigenVectors * eigenValuesSqrtInv * eigenVectors';

            % Apply Löwdin orthogonalization to original atomic orbitals
            % This transforms non-orthogonal atomic orbitals into orthogonal ones
            orthogonalizedOrbitals = symmetricOrthogMatrix * this.originalAtomicWavefunction;

            c = X' * orthogonalizedOrbitals;
        end

        function printAtomicStates(this)
            arguments
                this 
            end

        end
    end
end

