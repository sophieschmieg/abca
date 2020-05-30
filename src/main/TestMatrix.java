package main;

import fields.Element;
import fields.FreeModule;
import fields.Matrix;
import fields.MatrixAlgebra;
import fields.Polynomial;
import fields.Ring;
import fields.Vector;

public class TestMatrix<T extends Element> {
	public TestMatrix(Ring<T> ring, int dimension) {
		MatrixAlgebra<T> algebra = new MatrixAlgebra<>(ring, dimension);
		System.out.println(algebra.one());
		Matrix<T> random = algebra.getRandomElement();
		System.out.println(random);
		MatrixAlgebra<T>.GaussianEliminationResult decomp = algebra.gaussianElimination(random);
		System.out.println("R^-1.A.C^-1");
		Matrix<T> rac = algebra.multiply(decomp.getRowOperationsInverse(), random, decomp.getColOperationsInverse());
		System.out.println(rac);
		System.out.println(rac.equals(decomp.getDiagonalMatrix()));
		System.out.println("R.D.C");
		Matrix<T> rdc = algebra.multiply(decomp.getRowOperations(), decomp.getDiagonalMatrix(),
				decomp.getColOperations());
		System.out.println(rdc);
		System.out.println(rdc.equals(random));
		System.out.println("R.R^-1");
		Matrix<T> rr = algebra.multiply(decomp.getRowOperations(), decomp.getRowOperationsInverse());
		System.out.println(rr);
		System.out.println(rr.equals(algebra.one()));
		System.out.println("C^-1.C");
		Matrix<T> cc = algebra.multiply(decomp.getColOperationsInverse(), decomp.getColOperations());
		System.out.println(cc);
		System.out.println(cc.equals(algebra.one()));
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
		try {
			Polynomial<T> ch = algebra.characteristicPolynomial(random);
			System.out.println(ch);
			System.out.println(ch.getRing().factorization(ch));
		} catch (Exception e) {
		}
		FreeModule<T> free = new FreeModule<T>(ring, dimension);
		Vector<T> x = free.getRandomElement();
		System.out.println("Ker(A)");
		System.out.println(algebra.kernelBasis(random));
		System.out.println("x");
		System.out.println(x);
		System.out.println("A.x");
		Vector<T> b = algebra.multiply(random, x);
		System.out.println(b);
		System.out.println("A.x = b");
		System.out.println(algebra.solve(random, b));
	}
}
