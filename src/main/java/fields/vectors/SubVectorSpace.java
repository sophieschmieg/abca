package fields.vectors;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ring;
import fields.interfaces.VectorSpace;

public class SubVectorSpace<T extends Element<T>, S extends Element<S>> extends AbstractModule<T, S>
		implements VectorSpace<T, S> {
	private VectorSpace<T, S> space;
	private List<S> generators;
	private FiniteVectorSpace<T> asFiniteVectorSpace;

	public SubVectorSpace(VectorSpace<T, S> space, List<S> generators) {
		this.space = space;
		this.generators = new ArrayList<>();
		for (S s : generators) {
			this.generators.add(s);
			if (!space.isLinearIndependent(this.generators)) {
				this.generators.remove(this.generators.size() - 1);
			}
		}
		this.asFiniteVectorSpace = new FiniteVectorSpace<>(space.getField(), this.generators.size());
	}

	@Override
	public Exactness exactness() {
		return space.exactness();
	}

	@Override
	public Ring<T> getRing() {
		return space.getRing();
	}

	@Override
	public S zero() {
		return space.zero();
	}

	@Override
	public S add(S s1, S s2) {
		return space.add(s1, s2);
	}

	@Override
	public S negative(S s) {
		return space.negative(s);
	}

	@Override
	public S scalarMultiply(T t, S s) {
		return space.scalarMultiply(t, s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public boolean isLinearIndependent(List<S> s) {
		return space.isLinearIndependent(s);
	}

	@Override
	public List<List<T>> nonTrivialCombinations(List<S> s) {
		return space.nonTrivialCombinations(s);
	}

	public boolean contains(S s) {
		List<S> added = new ArrayList<>();
		added.addAll(generators);
		added.add(s);
		return !space.isLinearIndependent(added);
	}

	public SubVectorSpace<T, S> intersection(SubVectorSpace<T, S> other) {
		List<S> generators = new ArrayList<>();
		generators.addAll(this.generators);
		generators.addAll(other.generators);
		List<List<T>> coeffs = space.nonTrivialCombinations(generators);
		List<S> intersectionGenerators = new ArrayList<>();
		for (List<T> coeff : coeffs) {
			S generator = zero();
			for (int i = 0; i < this.generators.size(); i++) {
				generator = add(generator, scalarMultiply(coeff.get(i), this.generators.get(i)));
			}
			intersectionGenerators.add(generator);
		}
		return new SubVectorSpace<>(space, intersectionGenerators);
	}

	public SubVectorSpace<T, S> add(SubVectorSpace<T, S> other) {
		List<S> generators = new ArrayList<>();
		generators.addAll(this.generators);
		generators.addAll(other.generators);
		return new SubVectorSpace<>(space, generators);
	}

	@Override
	public boolean isGeneratingModule(List<S> s) {
		if (s.size() < dimension()) {
			return false;
		}
		List<S> independent = new ArrayList<>();
		for (S v : s) {
			if (!contains(v)) {
				return false;
			}
			independent.add(v);
			if (!space.isLinearIndependent(independent)) {
				independent.remove(independent.size() - 1);
			}
		}
		return independent.size() == dimension();
	}

	@Override
	public List<S> getModuleGenerators() {
		return this.generators;
	}

	@Override
	public S getRandomElement() {
		Vector<T> rng = asFiniteVectorSpace.getRandomElement();
		return fromVector(rng);
	}

	@Override
	public boolean isFinite() {
		return space.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return getField().getNumberOfElements().pow(generators.size());
	}

	@Override
	public Iterator<S> iterator() {
		return new Iterator<S>() {
			private Iterator<Vector<T>> it = asFiniteVectorSpace.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public S next() {
				return fromVector(it.next());
			}

		};
	}

	@Override
	public Field<T> getField() {
		return space.getField();
	}

	@Override
	public List<S> getBasis() {
		return generators;
	}

	@Override
	public Vector<T> asVector(S s) {
		List<S> columns = new ArrayList<>();
		columns.addAll(generators);
		columns.add(s);
		List<T> nonTrivial = space.nonTrivialCombination(columns);
		T last = nonTrivial.get(dimension());
		if (last.equals(space.getField().zero())) {
			throw new ArithmeticException("generators not linear independent!");
		}
		Vector<T> asVector = new Vector<>(nonTrivial.subList(0, dimension()));
		return asFiniteVectorSpace.scalarMultiply(space.getField().negative(space.getField().inverse(last)), asVector);
	}

	@Override
	public MatrixAlgebra<T> matrixAlgebra() {
		return asFiniteVectorSpace.matrixAlgebra();
	}

	@Override
	public int dimension() {
		return this.generators.size();
	}

	@Override
	public String toString() {
		StringBuilder build = new StringBuilder();
		build.append("<");
		boolean first = true;
		for (S generator : generators) {
			if (first) {
				first = false;
			} else {
				build.append(", ");
			}
			build.append(generator.toString());
		}
		build.append(">");
		return build.toString();
	}
}
