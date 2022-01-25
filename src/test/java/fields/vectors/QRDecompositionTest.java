package fields.vectors;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.InnerProductSpace.OrthogonalSimilarResult;
import fields.interfaces.InnerProductSpace.QRDecompositionResult;
import fields.interfaces.InnerProductSpace.SingularValueDecompositionResult;
import fields.interfaces.MathMap;

class QRDecompositionTest {

	private void printMatrix(Matrix<Real> t) {
		Reals r = Reals.r(10);
		Matrix<Real> rounded = Matrix.mapMatrix(new MathMap<Real, Real>() {
			@Override
			public Real evaluate(Real t) {
				return r.roundToFixedPoint(r.getEmbedding(t), 10);
			}
		}, t);
		System.out.println(rounded);
	}

	// @Test
	void test1() {
		Reals r = Reals.r(128);
		List<List<Real>> testVector = new ArrayList<>();
		List<Real> row1 = new ArrayList<>();
		row1.add(r.getInteger(12));
		row1.add(r.getInteger(-51));
		row1.add(r.getInteger(4));
		testVector.add(row1);
		List<Real> row2 = new ArrayList<>();
		row2.add(r.getInteger(6));
		row2.add(r.getInteger(167));
		row2.add(r.getInteger(-68));
		testVector.add(row2);
		List<Real> row3 = new ArrayList<>();
		row3.add(r.getInteger(-4));
		row3.add(r.getInteger(24));
		row3.add(r.getInteger(-41));
		testVector.add(row3);
		InnerProductSpace<Real, Vector<Real>> space = new FiniteRealVectorSpace(r, 3);
		MatrixAlgebra<Real> algebra = space.matrixAlgebra();
		Matrix<Real> t = new Matrix<>(testVector);
		printMatrix(t);
		printMatrix(algebra.inverse(t));
		printMatrix(algebra.multiply(algebra.inverse(t), t));
		Matrix<Real> hermitean = algebra.multiply(space.conjugateTranspose(t), t);
		printMatrix(hermitean);
		OrthogonalSimilarResult<Real> hessenberg = space.hessenbergForm(hermitean);
		QRDecompositionResult<Real> qrHessen = space.qrDecomposition(hessenberg.getOrthogonallySimilarMatrix(), true);
		printMatrix(
				algebra.multiply(space.conjugateTranspose(qrHessen.getUnitaryMatrix()), qrHessen.getUnitaryMatrix()));
		printMatrix(hessenberg.getOrthogonallySimilarMatrix());
		printMatrix(algebra.multiply(qrHessen.getUnitaryMatrix(), qrHessen.getUpperTriangularMatrix()));
		OrthogonalSimilarResult<Real> hermiteanSchurr = space.schurrForm(hermitean);
		printMatrix(hermiteanSchurr.getOrthogonallySimilarMatrix());
		QRDecompositionResult<Real> qr = space.qrDecomposition(hermitean);
		printMatrix(algebra.multiply(space.conjugateTranspose(qr.getUnitaryMatrix()), hermitean));
		printMatrix(qr.getUpperTriangularMatrix());
		OrthogonalSimilarResult<Real> schurr = space.schurrForm(t);
		printMatrix(schurr.getUnitaryMatrix());
		printMatrix(schurr.getOrthogonallySimilarMatrix());
		printMatrix(algebra.multiply(space.conjugateTranspose(schurr.getUnitaryMatrix()), schurr.getUnitaryMatrix()));
		printMatrix(algebra.multiply(space.conjugateTranspose(schurr.getUnitaryMatrix()),
				schurr.getOrthogonallySimilarMatrix(), schurr.getUnitaryMatrix()));
		SingularValueDecompositionResult<Real> svd = space.singularValueDecomposition(t);
		printMatrix(algebra.multiply(space.conjugateTranspose(svd.getLeftUnitaryMatrix()), svd.getDiagonalMatrix(),
				svd.getRightUnitaryMatrix()));
	}

	// @Test
	void test2() {
		Reals r = Reals.r(128);
		List<List<Real>> testVector = new ArrayList<>();
		List<Real> row1 = new ArrayList<>();
		row1.add(r.getInteger(12));
		row1.add(r.getInteger(-51));
		row1.add(r.getInteger(4));
		testVector.add(row1);
		List<Real> row2 = new ArrayList<>();
		row2.add(r.getInteger(6));
		row2.add(r.getInteger(167));
		row2.add(r.getInteger(-68));
		testVector.add(row2);
		testVector.add(row1);
		InnerProductSpace<Real, Vector<Real>> space = new FiniteRealVectorSpace(r, 3);
		MatrixAlgebra<Real> algebra = space.matrixAlgebra();
		Matrix<Real> t = new Matrix<>(testVector);
		printMatrix(t);
		Matrix<Real> hermitean = algebra.multiply(space.conjugateTranspose(t), t);
		printMatrix(hermitean);
		QRDecompositionResult<Real> qr = space.qrDecomposition(hermitean);
		printMatrix(algebra.multiply(space.conjugateTranspose(qr.getUnitaryMatrix()), hermitean));
		printMatrix(qr.getUpperTriangularMatrix());
		OrthogonalSimilarResult<Real> schurr = space.schurrForm(t);
		printMatrix(schurr.getUnitaryMatrix());
		printMatrix(schurr.getOrthogonallySimilarMatrix());
		printMatrix(algebra.multiply(space.conjugateTranspose(schurr.getUnitaryMatrix()), schurr.getUnitaryMatrix()));
		printMatrix(algebra.multiply(space.conjugateTranspose(schurr.getUnitaryMatrix()),
				schurr.getOrthogonallySimilarMatrix(), schurr.getUnitaryMatrix()));
		SingularValueDecompositionResult<Real> svd = space.singularValueDecomposition(t);
		printMatrix(algebra.multiply(space.conjugateTranspose(svd.getLeftUnitaryMatrix()), svd.getDiagonalMatrix(),
				svd.getRightUnitaryMatrix()));
		System.out.println(svd.getRank());
		System.out.println(svd.getDiagonalMatrix());
	}

	// @Test
	void test3() {
		Reals r = Reals.r(128);
		List<Real> diag = new ArrayList<>();
		diag.add(r.getInteger(1));
		diag.add(r.getInteger(3));
		diag.add(r.getInteger(2));
		InnerProductSpace<Real, Vector<Real>> space = new FiniteRealVectorSpace(r, 3);
		MatrixAlgebra<Real> algebra = space.matrixAlgebra();
		Matrix<Real> t = Matrix.fromDiagonal(diag, r.zero());
		printMatrix(t);
		SingularValueDecompositionResult<Real> svd = space.singularValueDecomposition(t);
		printMatrix(algebra.multiply(space.conjugateTranspose(svd.getLeftUnitaryMatrix()), svd.getDiagonalMatrix(),
				svd.getRightUnitaryMatrix()));

	}

	// @Test
	void test4() {
		Reals r = Reals.r(128);
		List<List<Real>> testVector = new ArrayList<>();
		List<Real> row1 = new ArrayList<>();
		row1.add(r.getInteger(0));
		row1.add(r.getInteger(-1));
		testVector.add(row1);
		List<Real> row2 = new ArrayList<>();
		row2.add(r.getInteger(1));
		row2.add(r.getInteger(0));
		testVector.add(row2);
		InnerProductSpace<Real, Vector<Real>> space = new FiniteRealVectorSpace(r, 2);
		Matrix<Real> t = new Matrix<>(testVector);
		MatrixModule<Real> algebra = t.getModule(r);
		printMatrix(t);
		SingularValueDecompositionResult<Real> svd = space.singularValueDecomposition(t);
		printMatrix(algebra.multiply(space.conjugateTranspose(svd.getLeftUnitaryMatrix()),
				algebra.multiply(svd.getDiagonalMatrix(), svd.getRightUnitaryMatrix())));
		printMatrix(svd.getDiagonalMatrix());
		printMatrix(svd.getLeftUnitaryMatrix());
		printMatrix(svd.getRightUnitaryMatrix());
	}

	// @Test
	void test5() {
		Reals r = Reals.r(128);
		List<List<Real>> testVector = new ArrayList<>();
		List<Real> row1 = new ArrayList<>();
		row1.add(r.getInteger(0));
		row1.add(r.getInteger(-1));
		row1.add(r.getInteger(-1));
		testVector.add(row1);
		List<Real> row2 = new ArrayList<>();
		row2.add(r.getInteger(1));
		row2.add(r.getInteger(0));
		row2.add(r.getInteger(1));
		testVector.add(row2);
		InnerProductSpace<Real, Vector<Real>> space = new FiniteRealVectorSpace(r, 2);
		Matrix<Real> t = new Matrix<>(testVector);
		MatrixModule<Real> module = t.getModule(r);
		printMatrix(t);
		SingularValueDecompositionResult<Real> svd = space.singularValueDecomposition(t);
		printMatrix(module.multiply(space.conjugateTranspose(svd.getLeftUnitaryMatrix()),
				module.multiply(svd.getDiagonalMatrix(), svd.getRightUnitaryMatrix())));
		printMatrix(svd.getDiagonalMatrix());
		printMatrix(svd.getLeftUnitaryMatrix());
		printMatrix(module.codomainAlgebra().multiply(space.conjugateTranspose(svd.getLeftUnitaryMatrix()),
				svd.getLeftUnitaryMatrix()));
		printMatrix(svd.getRightUnitaryMatrix());
		printMatrix(module.domainAlgebra().multiply(space.conjugateTranspose(svd.getRightUnitaryMatrix()),
				svd.getRightUnitaryMatrix()));
	}

	// @Test
	void test6() {
		Reals r = Reals.r(128);
		List<List<Real>> testVector = new ArrayList<>();
		List<Real> row1 = new ArrayList<>();
		row1.add(r.getInteger(0));
		row1.add(r.getInteger(-1));
		testVector.add(row1);
		List<Real> row2 = new ArrayList<>();
		row2.add(r.getInteger(1));
		row2.add(r.getInteger(0));
		testVector.add(row2);
		List<Real> row3 = new ArrayList<>();
		row3.add(r.getInteger(2));
		row3.add(r.getInteger(3));
		testVector.add(row3);
		InnerProductSpace<Real, Vector<Real>> space = new FiniteRealVectorSpace(r, 2);
		Matrix<Real> t = new Matrix<>(testVector);
		MatrixModule<Real> module = t.getModule(r);
		printMatrix(t);
		SingularValueDecompositionResult<Real> svd = space.singularValueDecomposition(t);
		printMatrix(module.multiply(space.conjugateTranspose(svd.getLeftUnitaryMatrix()),
				module.multiply(svd.getDiagonalMatrix(), svd.getRightUnitaryMatrix())));
		printMatrix(svd.getDiagonalMatrix());
		printMatrix(svd.getLeftUnitaryMatrix());
		printMatrix(module.codomainAlgebra().multiply(space.conjugateTranspose(svd.getLeftUnitaryMatrix()),
				svd.getLeftUnitaryMatrix()));
		printMatrix(svd.getRightUnitaryMatrix());
		printMatrix(module.domainAlgebra().multiply(space.conjugateTranspose(svd.getRightUnitaryMatrix()),
				svd.getRightUnitaryMatrix()));
	}

	@Test
	void test7() {
		Reals r = Reals.r(128);
		List<List<Real>> testVector = new ArrayList<>();
		List<Real> row1 = new ArrayList<>();
		row1.add(r.getInteger(1));
		row1.add(r.getInteger(-1));
		testVector.add(row1);
		List<Real> row2 = new ArrayList<>();
		row2.add(r.getInteger(-1));
		row2.add(r.getInteger(1));
		testVector.add(row2);
		InnerProductSpace<Real, Vector<Real>> space = new FiniteRealVectorSpace(r, 2);
		Matrix<Real> t = new Matrix<>(testVector);
		OrthogonalSimilarResult<Real> schurr = space.schurrForm(t);
		printMatrix(t);
		printMatrix(schurr.getOrthogonallySimilarMatrix());
		printMatrix(space.matrixAlgebra().multiply(space.conjugateTranspose(schurr.getUnitaryMatrix()),
				schurr.getOrthogonallySimilarMatrix(), schurr.getUnitaryMatrix()));
	}
}
