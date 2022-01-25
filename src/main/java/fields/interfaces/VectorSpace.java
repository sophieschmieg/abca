package fields.interfaces;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.MatrixModule;
import fields.vectors.SubVectorSpace;
import fields.vectors.Vector;

public interface VectorSpace<T extends Element<T>, S extends Element<S>> extends Module<T, S> {
	public Field<T> getField();

	public S getUnitVector(int index);

	public List<S> getBasis();

	public default boolean isBasis(List<S> s) {
		return s.size() == dimension() && isLinearIndependent(s);
	}

	public int dimension();

	public Vector<T> asVector(S s);

	public MatrixAlgebra<T> matrixAlgebra();

	@Override
	default List<S> getModuleGenerators() {
		return getBasis();
	}

	@Override
	default List<List<T>> getModuleGeneratorRelations() {
		return Collections.emptyList();
	}

	public default List<S> linearIndependentSubSet(List<S> s) {
		SubVectorSpace<T, S> asSubVectorSpace = new SubVectorSpace<>(this, s);
		List<Vector<T>> basis = new ArrayList<>();
		for (S basisVector : asSubVectorSpace.getBasis()) {
			Vector<T> asVector = asVector(basisVector);
			if (asVector.dimension() != asSubVectorSpace.getEmbeddingDimension()) {
				List<T> vectorList = new ArrayList<>();
				for (int i = 0; i < asSubVectorSpace.getEmbeddingDimension() - asVector.dimension(); i++) {
					vectorList.add(getField().zero());
				}
				vectorList.addAll(asVector.asList());
				asVector = new Vector<>(vectorList);
			}
			basis.add(asVector);
		}
		List<S> result = new ArrayList<>();
		Matrix<T> baseChange = Matrix.fromColumns(basis);
		MatrixModule<T> matrixModule = baseChange.getModule(getField());
		for (int i = 0; i < basis.size(); i++) {
			for (S t : s) {
				Vector<T> asVector = asVector(t);
				if (asVector.dimension() != asSubVectorSpace.getEmbeddingDimension()) {
					List<T> vectorList = new ArrayList<>();
					for (int j = 0; j < asSubVectorSpace.getEmbeddingDimension() - asVector.dimension(); j++) {
						vectorList.add(getField().zero());
					}
					vectorList.addAll(asVector.asList());
					asVector = new Vector<>(vectorList);
				}
				Vector<T> inBasis = matrixModule.solve(baseChange, asVector);
				if (inBasis.get(i + 1).equals(getField().zero())) {
					continue;
				}
				if (!getField().exactness().equals(Exactness.EXACT) && this instanceof InnerProductSpace<?, ?>) {
					InnerProductSpace<T, S> innerProduceSpace = (InnerProductSpace<T, S>) this;
					Reals r = innerProduceSpace.getValueField().getReals();
					Real epsilon = r.getPowerOfTwo(-r.precision() / 2);
					if (innerProduceSpace.getValueField().value(inBasis.get(i + 1)).compareTo(epsilon) < 0) {
						continue;
					}
				}
				basis.set(i, asVector);
				result.add(t);
				baseChange = Matrix.fromColumns(basis);
				break;
			}
		}
		return result;
	}
}
