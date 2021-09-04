package fields.interfaces;

import java.util.List;

import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

public interface VectorSpace<T extends Element<T>, S extends Element<S>> extends Module<T, S> {
	public Field<T> getField();
	public List<S> getBasis();
	public default boolean isBasis(List<S> s) {
		return s.size() == dimension() && isLinearIndependent(s);
	}
	public int dimension();
	public Vector<T> asVector(S s);
	public MatrixAlgebra<T> matrixAlgebra();
}
