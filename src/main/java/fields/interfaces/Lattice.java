package fields.interfaces;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Integers.IntE;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;

public interface Lattice<T extends Element<T>, B extends Element<B>, V extends Element<V>> extends Module<IntE, T> {
	public RealInnerProductSpace<B, V> getVectorSpace();

	public V embedding(T t);

	public int rank();

	public Matrix<B> generatorsAsMatrix();

	default public Real determinant() {
		Matrix<B> generators = generatorsAsMatrix();
		MatrixAlgebra<B> algebra = getVectorSpace().matrixAlgebra();
		if (generators.rows() == generators.columns()) {
			return getVectorSpace().asReal(algebra.determinant(generators));
		}
		Matrix<B> square = algebra.multiply(generators, algebra.transpose(generators));
		B squaredDeterminant = algebra.determinant(square);
		Reals r = getVectorSpace().getValueField().getReals();
		return r.positiveSqrt(getVectorSpace().asReal(squaredDeterminant));
	}
	
	default public Real rootHermiteFactor(T vector) {
		Reals r = getVectorSpace().getValueField().getReals();
		Real value = getVectorSpace().valueNorm(embedding(vector));
		Real root = r.positiveRoot(value, rank());
		return r.divide(root, determinant());
	}
}
