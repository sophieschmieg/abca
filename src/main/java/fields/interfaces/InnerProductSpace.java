package fields.interfaces;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals.Real;
import fields.vectors.DualVectorSpace;
import fields.vectors.Matrix;

public interface InnerProductSpace<T extends Element<T>, S extends Element<S>> extends NormedVectorSpace<T, S> {
	public ComplexNumber asComplexNumber(T t);

	public T fromComplexNumber(ComplexNumber t);

	public T fromReal(Real r);

	// s1*.s2
	public T innerProduct(S s1, S s2);

	@Override
	public Real valueNorm(S s);

	public List<S> gramSchmidt(List<S> s);

	public List<S> normedGramSchmidt(List<S> s);

	public List<S> extendToOrthonormalBasis(List<S> s);

	public T conjugateScalar(T s);

	public S conjugateVector(S t);

	public SesquilinearForm<T, S> getSesquilinearForm(Matrix<T> t);

	public static class SesquilinearFormDiagonalizationResult<T extends Element<T>, S extends Element<S>> {
		private List<S> basis;
		private int positive;
		private int negative;
		private int zero;
		private Matrix<T> diagonalMatrix;

		public SesquilinearFormDiagonalizationResult(ValueField<T> field, List<S> basis, int positive, int negative, int zero) {
			this.basis = basis;
			if (positive + negative + zero != basis.size()) {
				throw new ArithmeticException("dimensions mismatch");
			}
			List<T> diagonal = new ArrayList<>();
			for (int i = 0; i < positive; i++) {
				diagonal.add(field.one());
			}
			for (int i = 0; i < negative; i++) {
				diagonal.add(field.getInteger(-1));
			}
			for (int i = 0; i < zero; i++) {
				diagonal.add(field.zero());
			}
			this.diagonalMatrix = Matrix.fromDiagonal(diagonal, field.zero());
			this.positive = positive;
			this.negative = negative;
			this.zero = zero;
		}

		public List<S> getBasis() {
			return basis;
		}

		public Matrix<T> getDiagonalMatrix() {
			return diagonalMatrix;
		}

		public int getPositive() {
			return positive;
		}

		public int getNegative() {
			return negative;
		}

		public int getZero() {
			return zero;
		}
	}

	public SesquilinearFormDiagonalizationResult<T, S> diagonalizeSesquilinearForm(SesquilinearForm<T, S> form);

	public SesquilinearFormDiagonalizationResult<T, S> diagonalizeSesquilinearForm(Matrix<T> matrix);

	public DualVectorSpace<T, S> getDual();

	public InnerProductSpace<T, S> withDimension(int dimension);

	public Matrix<T> conjugateTranspose(Matrix<T> t);

	public class QRDecompositionResult<T extends Element<T>> {
		private Matrix<T> unitaryMatrix;
		private Matrix<T> upperTriangularMatrix;

		public QRDecompositionResult(Matrix<T> unitaryMatrix, Matrix<T> upperTriangularMatrix) {
			this.unitaryMatrix = unitaryMatrix;
			this.upperTriangularMatrix = upperTriangularMatrix;
		}

		public Matrix<T> getUnitaryMatrix() {
			return unitaryMatrix;
		}

		public Matrix<T> getUpperTriangularMatrix() {
			return upperTriangularMatrix;
		}

	}

	public QRDecompositionResult<T> qrDecomposition(Matrix<T> t);

	public QRDecompositionResult<T> qrDecomposition(Matrix<T> t, boolean hessenberg);

	public class OrthogonalSimilarResult<T extends Element<T>> {
		private Matrix<T> unitaryMatrix;
		private Matrix<T> orthogonallySimilarMatrix;

		public OrthogonalSimilarResult(Matrix<T> unitaryMatrix, Matrix<T> orthogonallySimilarMatrix) {
			this.unitaryMatrix = unitaryMatrix;
			this.orthogonallySimilarMatrix = orthogonallySimilarMatrix;
		}

		public Matrix<T> getUnitaryMatrix() {
			return unitaryMatrix;
		}

		public Matrix<T> getOrthogonallySimilarMatrix() {
			return orthogonallySimilarMatrix;
		}

	}

	public OrthogonalSimilarResult<T> schurForm(Matrix<T> t);

	public OrthogonalSimilarResult<T> hessenbergForm(Matrix<T> t);

	public class SingularValueDecompositionResult<T extends Element<T>> {
		private int rank;
		private Matrix<T> leftUnitaryMatrix;
		private Matrix<T> diagonalMatrix;
		private Matrix<T> rightUnitaryMatrix;

		public SingularValueDecompositionResult(int rank, Matrix<T> leftUnitaryMatrix, Matrix<T> diagonalMatrix,
				Matrix<T> rightUnitaryMatrix) {
			this.rank = rank;
			this.leftUnitaryMatrix = leftUnitaryMatrix;
			this.diagonalMatrix = diagonalMatrix;
			this.rightUnitaryMatrix = rightUnitaryMatrix;
		}

		public int getRank() {
			return rank;
		}

		public Matrix<T> getLeftUnitaryMatrix() {
			return leftUnitaryMatrix;
		}

		public Matrix<T> getDiagonalMatrix() {
			return diagonalMatrix;
		}

		public Matrix<T> getRightUnitaryMatrix() {
			return rightUnitaryMatrix;
		}
	}

	public SingularValueDecompositionResult<T> singularValueDecomposition(Matrix<T> t);

	public Matrix<T> pseudoInverse(Matrix<T> t);

	public Optional<Real> conditionNumber(Matrix<T> t);

}
