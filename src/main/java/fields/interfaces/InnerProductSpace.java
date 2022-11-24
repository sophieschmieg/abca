package fields.interfaces;

import java.util.List;
import java.util.Optional;

import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals.Real;
import fields.vectors.DualVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public interface InnerProductSpace<T extends Element<T>, S extends Element<S>> extends NormedVectorSpace<T, S> {
	public ComplexNumber asComplexNumber(T t);

	public T fromComplexNumber(ComplexNumber t);

	public T fromReal(Real r);

	public T innerProduct(S s1, S s2);

	@Override
	public Real valueNorm(S s);

	public List<S> gramSchmidt(List<S> s);

	public List<S> normedGramSchmidt(List<S> s);

	public List<S> extendToOrthonormalBasis(List<S> s);

	public T conjugate(T s);

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

	public OrthogonalSimilarResult<T> schurrForm(Matrix<T> t);

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
	
	public boolean isSubModuleMember(MatrixModule<T> module, Matrix<T> m, Vector<T> b);

	public Vector<T> asSubModuleMember(MatrixModule<T> module, Matrix<T> m, Vector<T> b);

	public List<Vector<T>> syzygyProblem(MatrixModule<T> module, Matrix<T> m);

	public List<Vector<T>> simplifySubModuleGenerators(MatrixModule<T> module, Matrix<T> m);
}
