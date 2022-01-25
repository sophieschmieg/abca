package fields.vectors;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;

class LDUPTest {

	@Test
	void testIntegersSquareMatrices() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		MatrixAlgebra<IntE> alg = new FreeModule<>(z, 4).matrixAlgebra();
		MatrixAlgebra<Fraction> m = new FiniteVectorSpace<>(q, 4).matrixAlgebra();
		for (int tc = 0; tc < 100; tc++) {
			Matrix<IntE> test = alg.getRandomElement();
			MatrixModule<IntE>.LDUPResult r = alg.ldup(test);
			Matrix<Fraction> lower = Matrix.mapMatrix(q.getEmbeddingMap(), r.getLowerTriangle());
			Matrix<Fraction> inverseDiagonal = Matrix.mapMatrix(q.getEmbeddingMap(), r.getInverseDiagonal());
			Matrix<Fraction> upper = Matrix.mapMatrix(q.getEmbeddingMap(), r.getUpperTriangle());
			assertEquals(Matrix.mapMatrix(q.getEmbeddingMap(), test),
					m.permuteRows(r.getPermutation(), m.multiply(lower, m.inverse(inverseDiagonal), upper)));
		}
	}

	@Test
	void testIntegersLongMatrices() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		MatrixModule<IntE> module = new MatrixModule<>(z, 4, 6);
		MatrixModule<Fraction> m = new MatrixModule<>(q, 4, 6);
		MatrixAlgebra<Fraction> alg = new FiniteVectorSpace<>(q, 4).matrixAlgebra();
		for (int tc = 0; tc < 100; tc++) {
			Matrix<IntE> test = module.getRandomElement();
			MatrixModule<IntE>.LDUPResult r = module.ldup(test);
			Matrix<Fraction> lower = Matrix.mapMatrix(q.getEmbeddingMap(), r.getLowerTriangle());
			Matrix<Fraction> inverseDiagonal = Matrix.mapMatrix(q.getEmbeddingMap(), r.getInverseDiagonal());
			Matrix<Fraction> upper = Matrix.mapMatrix(q.getEmbeddingMap(), r.getUpperTriangle());
			assertEquals(Matrix.mapMatrix(q.getEmbeddingMap(), test),
					m.permuteRows(r.getPermutation(), m.multiply(alg.multiply(lower, alg.inverse(inverseDiagonal)), upper)));
		}
	}

	@Test
	void testIntegersTallMatrices() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		MatrixModule<IntE> module = new MatrixModule<>(z, 6, 4);
		MatrixModule<Fraction> m = new MatrixModule<>(q, 6, 4);
		MatrixAlgebra<Fraction> alg = new FiniteVectorSpace<>(q, 6).matrixAlgebra();
		for (int tc = 0; tc < 100; tc++) {
			Matrix<IntE> test = module.getRandomElement();
			MatrixModule<IntE>.LDUPResult r = module.ldup(test);
			Matrix<Fraction> lower = Matrix.mapMatrix(q.getEmbeddingMap(), r.getLowerTriangle());
			Matrix<Fraction> inverseDiagonal = Matrix.mapMatrix(q.getEmbeddingMap(), r.getInverseDiagonal());
			Matrix<Fraction> upper = Matrix.mapMatrix(q.getEmbeddingMap(), r.getUpperTriangle());
			assertEquals(Matrix.mapMatrix(q.getEmbeddingMap(), test),
					m.permuteRows(r.getPermutation(), m.multiply(alg.multiply(lower, alg.inverse(inverseDiagonal)), upper)));
		}
	}

}
