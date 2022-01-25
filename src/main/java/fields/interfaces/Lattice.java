package fields.interfaces;

import fields.integers.Integers.IntE;
import fields.vectors.Matrix;

public interface Lattice<T extends Element<T>, B extends Element<B>, V extends Element<V>> extends Module<IntE, T> {
	public InnerProductSpace<B, V> getVectorSpace();

	public V embedding(T t);
	
	public int rank();
	
	public Matrix<B> generatorsAsMatrix();
}
