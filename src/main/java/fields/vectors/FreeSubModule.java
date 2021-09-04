package fields.vectors;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Module;
import fields.interfaces.Ring;

public class FreeSubModule<T extends Element<T>, S extends Element<S>> extends AbstractModule<T, S> {
	private Module<T, S> module;
	private FreeModule<T> freeModule;
	private List<S> generators;
	private List<Vector<T>> basisAsVectors;
	private Matrix<T> baseChange;
	private FreeModule<T> asFreeModule;

	public FreeSubModule(Module<T, S> module, List<S> generators) {
		if (!module.isFree()) {
			throw new ArithmeticException("module is not free");
		}
		this.module = module;
		List<Vector<T>> asVectors = new ArrayList<>();
		int embeddingDimension = 0;
		boolean differentEmbeddingDimensions = false;
		for (S s : generators) {
			Vector<T> vector = module.asVector(s);
			asVectors.add(vector);
			if (embeddingDimension != 0 && embeddingDimension != vector.dimension()) {
				differentEmbeddingDimensions = true;
			}
			if (embeddingDimension < vector.dimension()) {
				embeddingDimension = vector.dimension();
			}
		}
		this.freeModule = new FreeModule<>(module.getRing(), embeddingDimension);
		if (differentEmbeddingDimensions) {
			asVectors.clear();
			for (S s : generators) {
				Vector<T> vector = module.asVector(s);
				List<T> asList = new ArrayList<>();
				asList.addAll(vector.asList());
				while (asList.size() < embeddingDimension) {
					asList.add(module.getRing().zero());
				}
				asVectors.add(new Vector<>(asList));
			}
		}
		Matrix<T> asMatrix = Matrix.fromColumns(asVectors);
		this.basisAsVectors = asMatrix.getModule(module.getRing()).imageBasis(asMatrix);
		this.baseChange = Matrix.fromColumns(basisAsVectors);
		this.generators = new ArrayList<>();
		for (Vector<T> basisVector : basisAsVectors) {
			this.generators.add(module.fromVector(basisVector));
		}
		this.asFreeModule = new FreeModule<>(module.getRing(), this.generators.size());
	}

	@Override
	public Exactness exactness() {
		return module.exactness();
	}

	@Override
	public Ring<T> getRing() {
		return module.getRing();
	}

	@Override
	public S zero() {
		return module.zero();
	}

	@Override
	public S add(S s1, S s2) {
		return module.add(s1, s2);
	}

	@Override
	public S negative(S s) {
		return module.negative(s);
	}

	@Override
	public S scalarMultiply(T t, S s) {
		return module.scalarMultiply(t, s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public boolean isLinearIndependent(List<S> s) {
		return module.isLinearIndependent(s);
	}

	@Override
	public List<List<T>> nonTrivialCombinations(List<S> s) {
		return module.nonTrivialCombinations(s);
	}

	public boolean contains(S s) {
		return freeModule.matrixAlgebra().hasSolution(baseChange, module.asVector(s));
	}

	public FreeSubModule<T, S> add(FreeSubModule<T, S> other) {
		List<S> generators = new ArrayList<>();
		generators.addAll(this.generators);
		generators.addAll(other.generators);
		return new FreeSubModule<>(module, generators);
	}

	@Override
	public boolean isGeneratingModule(List<S> s) {
		if (s.size() < rank()) {
			return false;
		}
		List<Vector<T>> asVectors = new ArrayList<>();
		for (S v : s) {
			if (!contains(v)) {
				return false;
			}
			asVectors.add(asVector(v));
		}
		return asFreeModule.isGeneratingModule(asVectors);
	}

	@Override
	public List<S> getModuleGenerators() {
		return this.generators;
	}

	@Override
	public S getRandomElement() {
		Vector<T> rng = asFreeModule.getRandomElement();
		return fromVector(rng);
	}

	@Override
	public boolean isFinite() {
		return module.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return getRing().getNumberOfElements().pow(generators.size());
	}

	@Override
	public Iterator<S> iterator() {
		return new Iterator<S>() {
			private Iterator<Vector<T>> it = asFreeModule.iterator();

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

	public List<S> getBasis() {
		return generators;
	}

	@Override
	public Vector<T> asVector(S s) {
		return freeModule.matrixAlgebra().solve(baseChange, module.asVector(s));
	}

	public MatrixAlgebra<T> matrixAlgebra() {
		return asFreeModule.matrixAlgebra();
	}

	public int rank() {
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
