package fields.vectors;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals.Real;
import fields.interfaces.Element;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.ValueField;

public class SubInnerProductSpace<T extends Element<T>, S extends Element<S>> extends AbstractInnerProductSpace<T, S> {
	private SubVectorSpace<T, S> asSubVectorSpace;
	private InnerProductSpace<T, S> innerProductSpace;

	public SubInnerProductSpace(InnerProductSpace<T, S> innerProductSpace, SubVectorSpace<T, S> asSubVectorSpace) {
		this.innerProductSpace = innerProductSpace;
		this.asSubVectorSpace = asSubVectorSpace;
	}

	public SubInnerProductSpace(InnerProductSpace<T, S> innerProductSpace, List<S> generators) {
		this(innerProductSpace, new SubVectorSpace<>(innerProductSpace, generators));
	}

	@Override
	public T innerProduct(S s1, S s2) {
		return innerProductSpace.innerProduct(s1, s2);
	}

	@Override
	public T fromComplexNumber(ComplexNumber t) {
		return innerProductSpace.fromComplexNumber(t);
	}

	@Override
	public ComplexNumber asComplexNumber(T t) {
		return innerProductSpace.asComplexNumber(t);
	}

	@Override
	public T fromReal(Real r) {
		return innerProductSpace.fromReal(r);
	}

	@Override
	public T conjugate(T s) {
		return innerProductSpace.conjugate(s);
	}
	
	@Override
	public InnerProductSpace<T, S> withDimension(int dimension) {
		return innerProductSpace.withDimension(dimension);
	}

	@Override
	public ValueField<T> getValueField() {
		return innerProductSpace.getValueField();
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
	public List<List<T>> nonTrivialCombinations(List<S> s) {
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

	public SubInnerProductSpace<T, S> add(SubInnerProductSpace<T, S> other) {
		return new SubInnerProductSpace<>(innerProductSpace, asSubVectorSpace.add(other.asSubVectorSpace));
	}

	public SubInnerProductSpace<T, S> intersection(SubInnerProductSpace<T, S> other) {
		return new SubInnerProductSpace<>(innerProductSpace, asSubVectorSpace.intersection(other.asSubVectorSpace));
	}

	public SubInnerProductSpace<T, S> orthogonalComplement() {
		List<S> basis = new ArrayList<>();
		basis.addAll(gramSchmidt(getBasis()));
		SubVectorSpace<T, S> subVectorSpace = asSubVectorSpace;
		for (int i = 0; i < innerProductSpace.dimension(); i++) {
			if (!subVectorSpace.contains(innerProductSpace.getUnitVector(i + 1))) {
				basis.add(innerProductSpace.getUnitVector(i + 1));
				subVectorSpace = new SubVectorSpace<>(innerProductSpace, basis);
			}
		}
		basis = innerProductSpace.gramSchmidt(basis);
		return new SubInnerProductSpace<>(innerProductSpace, basis.subList(dimension(), innerProductSpace.dimension()));
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
