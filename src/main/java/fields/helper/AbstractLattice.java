package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Integers.IntegerIdeal;
import fields.interfaces.Element;
import fields.interfaces.Lattice;
import fields.interfaces.RealInnerProductSpace;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public abstract class AbstractLattice<R extends Element<R>, T extends Element<T>, S extends Element<S>>
		extends AbstractModule<IntE, R> implements Lattice<R, T, S> {
	private RealInnerProductSpace<T, S> embeddingSpace;
	private FreeModule<IntE> asFreeModule;
	private List<R> generators;
	private Matrix<T> generatorsAsMatrix;
	private MatrixModule<T> matrixModule;

	protected AbstractLattice(RealInnerProductSpace<T, S> embeddingSpace, List<R> generators, boolean isBasis) {
		Integers z = Integers.z();
		this.embeddingSpace = embeddingSpace;
		this.asFreeModule = new FreeModule<>(z, generators.size());
		this.generators = embeddingSpace.latticeReduction(generators, this, isBasis);
		List<Vector<T>> asReducedVectors = new ArrayList<>();
		for (R generator : this.generators) {
			asReducedVectors.add(embeddingSpace.asVector(embedding(generator)));
		}
		this.generatorsAsMatrix = Matrix.fromColumns(asReducedVectors);
		this.matrixModule = generatorsAsMatrix.getModule(embeddingSpace.getField());
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public IntegerIdeal annihilator() {
		return Integers.z().getZeroIdeal();
	}

	private List<Vector<IntE>> asIntegerList(List<R> s) {
		List<Vector<IntE>> result = new ArrayList<>();
		for (R t : s) {
			result.add(asVector(t));
		}
		return result;
	}

	@Override
	public boolean isLinearIndependent(List<R> s) {
		return asFreeModule.isLinearIndependent(asIntegerList(s));
	}

	@Override
	public boolean isGeneratingModule(List<R> s) {
		return asFreeModule.isGeneratingModule(asIntegerList(s));
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<R> s) {
		return asFreeModule.nonTrivialCombinations(asIntegerList(s));
	}

	@Override
	public List<R> getModuleGenerators() {
		return generators;
	}

	@Override
	public List<Vector<IntE>> getSyzygies() {
		return Collections.emptyList();
	}

	@Override
	public R getRandomElement() {
		return fromVector(asFreeModule.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<R> iterator() {
		throw new InfinityException();
	}

	@Override
	public RealInnerProductSpace<T, S> getVectorSpace() {
		return embeddingSpace;
	}

	@Override
	public int rank() {
		return getModuleGenerators().size();
	}

	@Override
	public Matrix<T> generatorsAsMatrix() {
		return generatorsAsMatrix;
	}

	@Override
	public Vector<IntE> asVector(R s) {
		Vector<T> embedding = embeddingSpace.asVector(embedding(s));
		Vector<T> inBasis = matrixModule.solve(generatorsAsMatrix, embedding);
		List<IntE> result = new ArrayList<>();
		for (T coeff : inBasis.asList()) {
			result.add(embeddingSpace.asReal(coeff).round());
		}
		return new Vector<>(result);
	}

}
