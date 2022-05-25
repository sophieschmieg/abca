package fields.vectors;

import java.util.Iterator;
import java.util.List;

import fields.floatingpoint.Reals.Real;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.RealInnerProductSpace;
import fields.interfaces.ValueField;

public class SubRealInnerProductSpace<T extends Element<T>, S extends Element<S>>
		extends AbstractRealInnerProductSpace<T, S> {
	private SubVectorSpace<T, S> asSubVectorSpace;
	private RealInnerProductSpace<T, S> innerProductSpace;

	public SubRealInnerProductSpace(RealInnerProductSpace<T, S> innerProductSpace,
			SubVectorSpace<T, S> asSubVectorSpace) {
		this.innerProductSpace = innerProductSpace;
		this.asSubVectorSpace = asSubVectorSpace;
	}

	public SubRealInnerProductSpace(RealInnerProductSpace<T, S> innerProductSpace, List<S> generators,
			boolean checkBasis) {
		this(innerProductSpace, new SubVectorSpace<>(innerProductSpace, generators, checkBasis));
	}

	public SubRealInnerProductSpace(RealInnerProductSpace<T, S> innerProductSpace, List<S> generators) {
		this(innerProductSpace, generators, true);
	}
	
	@Override
	public String toString() {
		return asSubVectorSpace.toString();
	}

	@Override
	public T innerProduct(S s1, S s2) {
		return innerProductSpace.innerProduct(s1, s2);
	}

	@Override
	public IntE round(T t) {
		return innerProductSpace.round(t);
	}

	@Override
	public T fromReal(Real r) {
		return innerProductSpace.fromReal(r);
	}

	@Override
	public Real asReal(T t) {
		return innerProductSpace.asReal(t);
	}

	@Override
	public ValueField<T> getValueField() {
		return innerProductSpace.getValueField();
	}

	@Override
	public RealInnerProductSpace<T, S> withDimension(int dimension) {
		return innerProductSpace.withDimension(dimension);
	}

	@Override
	public S getUnitVector(int index) {
		return asSubVectorSpace.getUnitVector(index);
	}

	@Override
	public List<S> getBasis() {
		return asSubVectorSpace.getBasis();
	}

	@Override
	public int dimension() {
		return asSubVectorSpace.dimension();
	}

	@Override
	public Vector<T> asVector(S s) {
		return asSubVectorSpace.asVector(s);
	}

	@Override
	public MatrixAlgebra<T> matrixAlgebra() {
		return asSubVectorSpace.matrixAlgebra();
	}

	@Override
	public S zero() {
		return innerProductSpace.zero();
	}

	@Override
	public S add(S s1, S s2) {
		return innerProductSpace.add(s1, s2);
	}

	@Override
	public S negative(S s) {
		return innerProductSpace.negative(s);
	}

	@Override
	public S scalarMultiply(T t, S s) {
		return innerProductSpace.scalarMultiply(t, s);
	}

	@Override
	public boolean isGeneratingModule(List<S> s) {
		return asSubVectorSpace.isLinearIndependent(s);
	}

	@Override
	public List<Vector<T>> nonTrivialCombinations(List<S> s) {
		return asSubVectorSpace.nonTrivialCombinations(s);
	}

	@Override
	public Exactness exactness() {
		return innerProductSpace.exactness();
	}

	@Override
	public S getRandomElement() {
		return asSubVectorSpace.getRandomElement();
	}

	@Override
	public Iterator<S> iterator() {
		return asSubVectorSpace.iterator();
	}

	public boolean contains(S s) {
		return asSubVectorSpace.contains(s);
	}

	public SubRealInnerProductSpace<T, S> add(SubRealInnerProductSpace<T, S> other) {
		return new SubRealInnerProductSpace<>(innerProductSpace, asSubVectorSpace.add(other.asSubVectorSpace));
	}

	public SubRealInnerProductSpace<T, S> intersection(SubRealInnerProductSpace<T, S> other) {
		return new SubRealInnerProductSpace<>(innerProductSpace, asSubVectorSpace.intersection(other.asSubVectorSpace));
	}

	public SubRealInnerProductSpace<T, S> orthogonalComplement() {
		List<S> basis = innerProductSpace.extendToOrthonormalBasis(getBasis());
		return new SubRealInnerProductSpace<>(innerProductSpace,
				basis.subList(dimension(), innerProductSpace.dimension()), false);
	}

	public S project(S s) {
		List<S> basis = gramSchmidt(getBasis());
		S result = zero();
		for (S basisVector : basis) {
			result = add(scalarMultiply(
					getField().divide(innerProduct(s, basisVector), innerProduct(basisVector, basisVector)),
					basisVector), result);
		}
		return result;
	}
}
