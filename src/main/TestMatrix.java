package main;

import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ring;
import fields.polynomials.UnivariatePolynomial;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public class TestMatrix<T extends Element<T>> {
	public TestMatrix(Ring<T> ring, int rows, int cols) {
		MatrixModule<T> module = new MatrixModule<>(ring, rows, cols);
		Matrix<T> random = module.getRandomElement();
		System.out.println(random);
		if (!ring.isEuclidean()) {
			return;
		}
		MatrixModule<T>.GaussianEliminationResult decomp = module.gaussianElimination(random);
		System.out.println("R^-1.A.C^-1");
		Matrix<T> rac = module.multiply(module.multiply(decomp.getRowOperationsInverse(), random),
				decomp.getColOperationsInverse());
		System.out.println(rac);
		System.out.println(rac.equals(decomp.getDiagonalMatrix()));
		System.out.println("R.D.C");
		Matrix<T> rdc = module.multiply(module.multiply(decomp.getRowOperations(), decomp.getDiagonalMatrix()),
				decomp.getColOperations());
		System.out.println(rdc);
		System.out.println(rdc.equals(random));
		System.out.println("R.R^-1");
		Matrix<T> rr = module.codomainAlgebra().multiply(decomp.getRowOperations(), decomp.getRowOperationsInverse());
		System.out.println(rr);
		System.out.println(rr.equals(module.codomainAlgebra().one()));
		System.out.println("C^-1.C");
		Matrix<T> cc = module.domainAlgebra().multiply(decomp.getColOperationsInverse(), decomp.getColOperations());
		System.out.println(cc);
		System.out.println(cc.equals(module.domainAlgebra().one()));
		if (rows == cols) {
			MatrixAlgebra<T> algebra = module.domainAlgebra();
			System.out.println("det(A) = " + algebra.determinant(random));
			if (algebra.isUnit(random)) {
				System.out.println("A^-1");
				Matrix<T> inv = algebra.inverse(random);
				System.out.println(inv);
				System.out.println("A^-1.A");
				System.out.println(algebra.multiply(inv, random));
				System.out.println("A.A^-1");
				System.out.println(algebra.multiply(random, inv));
			}
			if (ring.getUnivariatePolynomialRing().isEuclidean()) {
				UnivariatePolynomial<T> ch = algebra.characteristicPolynomial(random);
				System.out.println(ch);
				if (ring instanceof Field) {
					System.out.println(((Field<T>) ring).factorization(ch));
				}
			}
		}
		FreeModule<T> free = module.domain();
		Vector<T> x = free.getRandomElement();
		System.out.println("Ker(A)");
		System.out.println(module.kernelBasis(random));
		System.out.println("x");
		System.out.println(x);
		System.out.println("A.x");
		Vector<T> b = module.multiply(random, x);
		System.out.println(b);
		System.out.println("A.x = b");
		Vector<T> solution = module.solve(random, b);
		System.out.println(solution);
		System.out.println(module.multiply(random, solution));
	}
}
