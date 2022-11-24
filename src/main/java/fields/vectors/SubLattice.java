package fields.vectors;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractModule;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Lattice;
import fields.interfaces.RealInnerProductSpace;

public class SubLattice<R extends Element<R>, T extends Element<T>, S extends Element<S>>
		extends AbstractModule<IntE, R> implements Lattice<R, T, S> {
	private Lattice<R, T, S> lattice;
	private List<R> generators;
	private Matrix<IntE> sublatticeBaseInLattice;
	private MatrixModule<IntE> matrixModule;
	private FreeSubModule<IntE, R> asFreeSubModule;
	private Matrix<T> generatorsAsMatrix;

	public SubLattice(Lattice<R, T, S> lattice, List<R> generators) {
		this.lattice = lattice;
		this.generators = lattice.getVectorSpace().latticeReduction(generators, lattice);
		this.asFreeSubModule = new FreeSubModule<>(lattice, this.generators);
		List<Vector<T>> asVectors = new ArrayList<>();
		for (R generator : this.generators) {
			asVectors.add(lattice.getVectorSpace().asVector(lattice.embedding(generator)));
		}
		this.generatorsAsMatrix = Matrix.fromColumns(asVectors);
		List<Vector<IntE>> inLatticeBase = new ArrayList<>();
		for (R base : this.generators) {
			inLatticeBase.add(lattice.asVector(base));
		}
		this.sublatticeBaseInLattice = Matrix.fromColumns(inLatticeBase);
		this.matrixModule = sublatticeBaseInLattice.getModule(Integers.z());
	}

	public SubLattice(Lattice<R, T, S> lattice, FreeSubModule<IntE, R> asFreeSubModule) {
		this(lattice, asFreeSubModule.getBasis());
	}

	@Override
	public String toString() {
		return generators.toString();
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public R zero() {
		return lattice.zero();
	}

	@Override
	public R add(R s1, R s2) {
		return lattice.add(s1, s2);
	}

	@Override
	public R negative(R s) {
		return lattice.negative(s);
	}

	@Override
	public R scalarMultiply(IntE t, R s) {
		return lattice.scalarMultiply(t, s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public Ideal<IntE> annihilator() {
		return Integers.z().getZeroIdeal();
	}

	@Override
	public List<Vector<IntE>> getSyzygies() {
		return Collections.emptyList();
	}

	@Override
	public boolean isLinearIndependent(List<R> s) {
		return asFreeSubModule.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<R> s) {
		return asFreeSubModule.isGeneratingModule(s);
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<R> s) {
		return asFreeSubModule.nonTrivialCombinations(s);
	}

	@Override
	public List<R> getModuleGenerators() {
		return generators;
	}

	@Override
	public Vector<IntE> asVector(R s) {
		Vector<IntE> inLatticeBase = lattice.asVector(s);
		return matrixModule.solve(sublatticeBaseInLattice, inLatticeBase);
	}

	@Override
	public Exactness exactness() {
		return lattice.exactness();
	}

	@Override
	public R getRandomElement() {
		return asFreeSubModule.getRandomElement();
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
		return asFreeSubModule.iterator();
	}

	@Override
	public RealInnerProductSpace<T, S> getVectorSpace() {
		return lattice.getVectorSpace();
	}

	@Override
	public S embedding(R t) {
		return lattice.embedding(t);
	}

	@Override
	public int rank() {
		return asFreeSubModule.rank();
	}

	@Override
	public Matrix<T> generatorsAsMatrix() {
		return generatorsAsMatrix;
	}

	public FreeSubModule<IntE, R> asFreeSubModule() {
		return asFreeSubModule;
	}
}
