package fields.vectors;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.Ring;
import fields.interfaces.VectorSpace;

public class SubVectorSpace<T extends Element<T>, S extends Element<S>> extends AbstractModule<T, S>
		implements VectorSpace<T, S> {
	private VectorSpace<T, S> space;
	private Matrix<T> baseChange;
	private MatrixModule<T> matrixModule;
	private List<S> generators;
	private FiniteVectorSpace<T> asFiniteVectorSpace;
	private int embeddingDimension;

	public SubVectorSpace(VectorSpace<T, S> space, List<S> generators, boolean checkBasis) {
		this.space = space;
		if (generators.size() == 0) {
			this.generators = Collections.emptyList();
			return;
		}
		List<Vector<T>> asVectors = new ArrayList<>();
		embeddingDimension = 0;
		boolean differentEmbeddingDimensions = false;
		for (S s : generators) {
			Vector<T> vector = space.asVector(s);
			asVectors.add(vector);
			if (embeddingDimension != 0 && embeddingDimension != vector.dimension()) {
				System.err.println("Different embeddingDimension: " + embeddingDimension + " != " + vector.dimension());
				differentEmbeddingDimensions = true;
			}
			if (embeddingDimension < vector.dimension()) {
				embeddingDimension = vector.dimension();
//			System.err.println("embeddingDimension is " + embeddingDimension);
			}
		}
		if (differentEmbeddingDimensions) {
			asVectors.clear();
			for (S s : generators) {
				Vector<T> vector = space.asVector(s);
				List<T> asList = new ArrayList<>();
				asList.addAll(vector.asList());
				while (asList.size() < embeddingDimension) {
					asList.add(space.getField().zero());
				}
				asVectors.add(new Vector<>(asList));
			}
		}
		if (asVectors.isEmpty() || embeddingDimension == 0) {
			System.err.println("Empty vector list " + asVectors + " or embeddingDimension 0: " + embeddingDimension);
			this.generators = Collections.emptyList();
			this.asFiniteVectorSpace = new FiniteVectorSpace<>(space.getField(), 0);
			return;
		}
		Matrix<T> asMatrix = Matrix.fromColumns(asVectors);
		this.matrixModule = asMatrix.getModule(space.getField());
		List<Vector<T>> basisVectors = matrixModule.imageBasis(asMatrix);
		if (basisVectors.isEmpty()) {
			System.err.println("No image basis found for: " + asVectors);
			this.generators = Collections.emptyList();
			this.asFiniteVectorSpace = new FiniteVectorSpace<>(space.getField(), 0);
			return;
		}
		this.baseChange = Matrix.fromColumns(basisVectors);
		this.matrixModule = baseChange.getModule(space.getField());
		this.generators = new ArrayList<>();
		for (Vector<T> basisVector : basisVectors) {
			this.generators.add(space.fromVector(basisVector));
		}
		this.asFiniteVectorSpace = new FiniteVectorSpace<>(space.getField(), this.generators.size());
	}

	public SubVectorSpace(VectorSpace<T, S> space, List<S> generators) {
		this(space, generators, true);
	}

	public int getEmbeddingDimension() {
		return embeddingDimension;
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
	public Ideal<T> annihilator() {
		return space.getRing().getZeroIdeal();
	}

	@Override
	public boolean isLinearIndependent(List<S> s) {
		return space.isLinearIndependent(s);
	}

	@Override
	public List<Vector<T>> nonTrivialCombinations(List<S> s) {
		return space.nonTrivialCombinations(s);
	}

	public boolean contains(S s) {
		if (generators.isEmpty()) {
			return s.equals(space.zero());
		}
		return matrixModule.hasSolution(baseChange, space.asVector(s));
	}

	public SubVectorSpace<T, S> intersection(SubVectorSpace<T, S> other) {
		if (this.generators.isEmpty()) {
			return this;
		}
		if (other.generators.isEmpty()) {
			return other;
		}
		List<S> generators = new ArrayList<>();
		generators.addAll(this.generators);
		generators.addAll(other.generators);
		List<Vector<T>> coeffs = space.nonTrivialCombinations(generators);
		List<S> intersectionGenerators = new ArrayList<>();
		for (Vector<T> coeff : coeffs) {
			S generator = zero();
			for (int i = 0; i < this.generators.size(); i++) {
				generator = add(generator, scalarMultiply(coeff.get(i + 1), this.generators.get(i)));
			}
			intersectionGenerators.add(generator);
		}
		return new SubVectorSpace<>(space, intersectionGenerators);
	}

	public SubVectorSpace<T, S> add(SubVectorSpace<T, S> other) {
		if (this.generators.isEmpty()) {
			return other;
		}
		if (other.generators.isEmpty()) {
			return this;
		}
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
		if (generators.isEmpty()) {
			return true;
		}
		List<Vector<T>> asVectors = new ArrayList<>();
		for (S t : s) {
			if (!contains(t)) {
				return false;
			}
			asVectors.add(space.asVector(t));
		}
		Matrix<T> asMatrix = Matrix.fromColumns(asVectors);
		int rank = asMatrix.getModule(space.getField()).rank(asMatrix);
		return rank == dimension();
	}

	@Override
	public List<S> getModuleGenerators() {
		return this.generators;
	}

	@Override
	public S getRandomElement() {
		if (generators.isEmpty()) {
			return space.zero();
		}
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
		if (generators.isEmpty()) {
			return Collections.singletonList(space.zero()).iterator();
		}
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
	public S getUnitVector(int index) {
		return generators.get(index - 1);
	}

	@Override
	public List<S> getBasis() {
		return generators;
	}

	@Override
	public Vector<T> asVector(S s) {
		if (generators.isEmpty()) {
			return new Vector<>();
		}
		return matrixModule.solve(baseChange, space.asVector(s));
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
