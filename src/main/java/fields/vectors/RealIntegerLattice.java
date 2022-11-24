package fields.vectors;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractModule;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Integers.IntegerIdeal;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Lattice;
import fields.interfaces.RealInnerProductSpace;

public class RealIntegerLattice extends AbstractModule<IntE, Vector<IntE>>
		implements Lattice<Vector<IntE>, Real, Vector<Real>> {
	private FreeModule<IntE> freeModule;
	private FreeModule<IntE> rankModule;
	private FiniteRealVectorSpace space;
	private List<Vector<IntE>> basis;
	private Matrix<Real> baseChange;
	private Matrix<IntE> integerBaseChange;
	private MatrixModule<IntE> integerMatrixModule;
	private Matrix<Fraction> rationalBaseChange;
	private MatrixModule<Fraction> rationalMatrixModule;

	public RealIntegerLattice(Reals r, int dimension, List<Vector<IntE>> basis) {
		this(new FiniteRealVectorSpace(r, dimension), new FreeModule<>(Integers.z(), dimension), basis);
	}

	public RealIntegerLattice(FiniteRealVectorSpace space, FreeModule<IntE> freeModule, List<Vector<IntE>> basis) {
		this.space = space;
		this.freeModule = freeModule;
		this.rankModule = new FreeModule<>(Integers.z(), basis.size());
		this.basis = space.latticeReduction(basis, this);
		if (basis.isEmpty()) {
			return;
		}
		List<Vector<Real>> embedded = new ArrayList<>();
		for (Vector<IntE> basisVector : this.basis) {
			embedded.add(embedding(basisVector));
		}
		this.baseChange = Matrix.fromColumns(embedded);
		this.integerBaseChange = Matrix.fromColumns(this.basis);
		this.integerMatrixModule = integerBaseChange.getModule(Integers.z());
		this.rationalBaseChange = Matrix.mapMatrix(Rationals.q().getEmbeddingMap(), integerBaseChange);
		this.rationalMatrixModule = rationalBaseChange.getModule(Rationals.q());
	}

	@Override
	public String toString() {
		StringBuilder build = new StringBuilder();
		build.append("<");
		boolean first = true;
		for (Vector<IntE> generator : basis) {
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

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public Vector<IntE> zero() {
		return freeModule.zero();
	}

	@Override
	public Vector<IntE> add(Vector<IntE> s1, Vector<IntE> s2) {
		return freeModule.add(s1, s2);
	}

	@Override
	public Vector<IntE> negative(Vector<IntE> s) {
		return freeModule.negative(s);
	}

	@Override
	public Vector<IntE> scalarMultiply(IntE t, Vector<IntE> s) {
		return freeModule.scalarMultiply(t, s);
	}

	@Override
	public boolean isFree() {
		return true;
	}
	
	public FreeModule<IntE> getFreeModule() {
		return freeModule;
	}

	public boolean contains(Vector<IntE> s) {
		return integerMatrixModule.hasSolution(integerBaseChange, s);
	}

	public boolean contains(FreeSubModule<IntE, Vector<IntE>> other) {
		for (Vector<IntE> otherBasisElement : other.getBasis()) {
			if (!contains(otherBasisElement)) {
				return false;
			}
		}
		return true;
	}

	public boolean contains(RealIntegerLattice other) {
		for (Vector<IntE> otherBasisElement : other.getBasis()) {
			if (!contains(otherBasisElement)) {
				return false;
			}
		}
		return true;
	}

	public FreeSubModule<IntE, Vector<IntE>> asFreeSubModule() {
		return new FreeSubModule<>(freeModule, basis);
	}

	public RealIntegerLattice intersection(RealIntegerLattice other) {
		List<Vector<IntE>> generators = new ArrayList<>();
		generators.addAll(this.basis);
		generators.addAll(other.basis);
		List<Vector<IntE>> coeffs = Integers.z().syzygyProblem(generators);
		List<Vector<IntE>> intersectionGenerators = new ArrayList<>();
		for (Vector<IntE> coeff : coeffs) {
			Vector<IntE> generator = zero();
			for (int i = 0; i < this.basis.size(); i++) {
				generator = add(generator, scalarMultiply(coeff.get(i + 1), this.basis.get(i)));
			}
			intersectionGenerators.add(generator);
		}
		return new RealIntegerLattice(space, freeModule, intersectionGenerators);
	}

	public RealIntegerLattice add(RealIntegerLattice other) {
		List<Vector<IntE>> generators = new ArrayList<>();
		generators.addAll(this.basis);
		generators.addAll(other.basis);
		return new RealIntegerLattice(space, freeModule,
				Integers.z().simplifySubModuleGenerators(generators));
	}

	@Override
	public IntegerIdeal annihilator() {
		return Integers.z().getZeroIdeal();
	}

	@Override
	public boolean isLinearIndependent(List<Vector<IntE>> s) {
		return freeModule.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<Vector<IntE>> s) {
		List<Vector<IntE>> inRankModule = new ArrayList<>();
		for (Vector<IntE> t : s) {
			inRankModule.add(asVector(t));
		}
		return rankModule.isGeneratingModule(inRankModule);
	}

	public boolean isBasis(List<Vector<IntE>> s) {
		return isGeneratingModule(s) && isLinearIndependent(s);
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<Vector<IntE>> s) {
		return freeModule.nonTrivialCombinations(s);
	}

	@Override
	public List<Vector<IntE>> getModuleGenerators() {
		return basis;
	}

	public List<Vector<IntE>> getBasis() {
		return basis;
	}

	public Vector<IntE> closestLatticePoint(Vector<IntE> t) {
		Vector<Fraction> solved = solveIntVector(t);
		List<IntE> coeffs = new ArrayList<>();
		for (Fraction coeff : solved.asList()) {
			coeffs.add(coeff.round());
		}
		return fromVector(new Vector<>(coeffs));
	}

	private Vector<Fraction> solveIntVector(Vector<IntE> b) {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		List<Fraction> result = new ArrayList<>();
		MatrixModule<IntE>.SmithNormalFormResult gauss = integerMatrixModule.smithNormalForm(integerBaseChange);
		Vector<IntE> rhs = integerMatrixModule.codomainAlgebra().multiply(gauss.getRowOperationsInverse(), b);
		for (int i = 0; i < baseChange.columns(); i++) {
			if (i < baseChange.rows() && (!rhs.get(i + 1).equals(z.zero())
					|| !gauss.getDiagonalMatrix().entry(i + 1, i + 1).equals(z.zero()))) {
				Fraction r = q.divide(q.getInteger(rhs.get(i + 1)),
						q.getInteger(gauss.getDiagonalMatrix().entry(i + 1, i + 1)));
				result.add(r);
			} else {
				result.add(q.zero());
			}
		}
		return rationalMatrixModule.domainAlgebra()
				.multiply(Matrix.mapMatrix(q.getEmbeddingMap(), gauss.getColOperationsInverse()), new Vector<>(result));
	}

	@Override
	public List<Vector<IntE>> getSyzygies() {
		return Collections.emptyList();
	}

	@Override
	public Vector<IntE> asVector(Vector<IntE> s) {
		return integerMatrixModule.solve(integerBaseChange, s);
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public Vector<IntE> getRandomElement() {
		return fromVector(rankModule.getRandomElement());
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
	public Iterator<Vector<IntE>> iterator() {
		return new Iterator<>() {
			private Iterator<Vector<IntE>> it = rankModule.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Vector<IntE> next() {
				return fromVector(it.next());
			}
		};
	}

	@Override
	public FiniteRealVectorSpace getVectorSpace() {
		return space;
	}

	@Override
	public Vector<Real> embedding(Vector<IntE> t) {
		Reals r = space.getValueField();
		List<Real> result = new ArrayList<>();
		for (IntE coeff : t.asList()) {
			result.add(r.getInteger(coeff));
		}
		return new Vector<>(result);
	}

	@Override
	public int rank() {
		return basis.size();
	}

	@Override
	public Matrix<Real> generatorsAsMatrix() {
		return baseChange;
	}

	public Vector<IntE> residue(Vector<IntE> t) {
		Vector<Real> embedded = embedding(t);
		Vector<IntE> closest = space.closestLatticePoint(embedded, this);
		return freeModule.subtract(t, closest);
	}

	public Matrix<IntE> integerBaseChange() {
		return integerBaseChange;
	}
	
	public Matrix<Fraction> rationalBaseChange() {
		return rationalBaseChange;
	}

}
