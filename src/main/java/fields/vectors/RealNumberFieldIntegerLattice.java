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
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers;

public class RealNumberFieldIntegerLattice extends AbstractModule<IntE, Vector<NFE>>
		implements Lattice<Vector<NFE>, Real, Vector<Real>> {
	private FreeModule<NFE> free;
	private FreeModule<IntE> asFreeIntegers;
	private FreeModule<IntE> embeddedFreeIntegers;
	private FreeSubModule<IntE, Vector<IntE>> asEmbeddedFreeIntegersSubModule;
	private NumberFieldIntegers order;
	private FiniteRealVectorSpace space;
	private List<Vector<NFE>> basis;
	private Matrix<Real> baseChange;
	private Matrix<IntE> integerBaseChange;
	private MatrixModule<IntE> integerMatrixModule;
	private Matrix<Fraction> rationalBaseChange;
	private MatrixModule<Fraction> rationalMatrixModule;

	public RealNumberFieldIntegerLattice(Reals r, int dimension, NumberFieldIntegers order, List<Vector<NFE>> basis) {
		this.order = order;
		this.space = new FiniteRealVectorSpace(r, dimension * order.rank());
		this.free = new FreeModule<>(order, dimension);
		this.basis = space.latticeReduction(basis, this);
		if (this.basis.isEmpty()) {
			return;
		}
		this.asFreeIntegers = new FreeModule<>(Integers.z(), this.basis.size());
		this.embeddedFreeIntegers = new FreeModule<>(Integers.z(), dimension * order.rank());
		List<Vector<Real>> embedded = new ArrayList<>();
		List<Vector<IntE>> integerEmbedded = new ArrayList<>();
		for (Vector<NFE> basisVector : this.basis) {
			embedded.add(embedding(basisVector));
			integerEmbedded.add(asIntegerVector(basisVector));
		}
		this.asEmbeddedFreeIntegersSubModule = new FreeSubModule<>(embeddedFreeIntegers, integerEmbedded);
		this.baseChange = Matrix.fromColumns(embedded);
		this.integerBaseChange = Matrix.fromColumns(integerEmbedded);
		this.integerMatrixModule = integerBaseChange.getModule(Integers.z());
		this.rationalBaseChange = Matrix.mapMatrix(Rationals.q().getEmbeddingMap(), integerBaseChange);
		this.rationalMatrixModule = rationalBaseChange.getModule(Rationals.q());
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public Vector<NFE> zero() {
		return free.zero();
	}

	@Override
	public Vector<NFE> add(Vector<NFE> s1, Vector<NFE> s2) {
		return free.add(s1, s2);
	}

	@Override
	public Vector<NFE> negative(Vector<NFE> s) {
		return free.negative(s);
	}

	@Override
	public Vector<NFE> scalarMultiply(IntE t, Vector<NFE> s) {
		return free.scalarMultiply(order.getInteger(t), s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public IntegerIdeal annihilator() {
		return Integers.z().getZeroIdeal();
	}

	private Vector<IntE> asIntegerVector(Vector<NFE> t) {
		List<IntE> result = new ArrayList<>();
		for (NFE c : t.asList()) {
			result.addAll(order.asVector(c).asList());
		}
		return new Vector<>(result);
	}
	
	private List<Vector<IntE>> asIntegerVectorList(List<Vector<NFE>> t) {
		List<Vector<IntE>> result = new ArrayList<>();
		for (Vector<NFE> v : t) {
			result.add(asIntegerVector(v));
		}
		return result;
	}

//	private Vector<NFE> fromIntegerVector(Vector<IntE> t) {
//		List<NFE> result = new ArrayList<>();
//		List<IntE> asList = t.asList();
//		for (int i = 0; i < free.dimension(); i++) {
//			result.add(order.fromVector(new Vector<>(asList.subList(i * order.rank(), (i + 1) * order.rank()))));
//		}
//		return new Vector<>(result);
//	}
//	
//	private List<Vector<NFE>> fromIntegerVectorList(List<Vector<IntE>> t) {
//		List<Vector<NFE>> result = new ArrayList<>();
//		for (Vector<IntE> v : t) {
//			result.add(fromIntegerVector(v));
//		}
//		return result;
//	}

	@Override
	public boolean isLinearIndependent(List<Vector<NFE>> s) {
		return free.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<Vector<NFE>> s) {
		return asEmbeddedFreeIntegersSubModule.isGeneratingModule(asIntegerVectorList(s));
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<Vector<NFE>> s) {
		return asEmbeddedFreeIntegersSubModule.nonTrivialCombinations(asIntegerVectorList(s));
	}

	@Override
	public List<Vector<NFE>> getModuleGenerators() {
		return basis;
	}

	public Vector<NFE> closestLatticePoint(Vector<NFE> t) {
		Vector<Fraction> solved = solveIntVector(asIntegerVector(t));
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
	public Vector<IntE> asVector(Vector<NFE> s) {
		return integerMatrixModule.solve(integerBaseChange, asIntegerVector(s));
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public Vector<NFE> getRandomElement() {
		return fromVector(asFreeIntegers.getRandomElement());
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
	public Iterator<Vector<NFE>> iterator() {
		return new Iterator<Vector<NFE>>() {
			private Iterator<Vector<IntE>> it = asFreeIntegers.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Vector<NFE> next() {
				return fromVector(it.next());
			}
		};
	}

	@Override
	public RealInnerProductSpace<Real, Vector<Real>> getVectorSpace() {
		return space;
	}

	@Override
	public Vector<Real> embedding(Vector<NFE> t) {
		List<Real> result = new ArrayList<>();
		for (NFE coeff : t.asList()) {
			result.addAll(order.embedding(coeff).asList());
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

	public Vector<NFE> residue(Vector<NFE> t) {
		Vector<Real> embedded = embedding(t);
		Vector<NFE> closest = space.closestLatticePoint(embedded, this);
		return free.subtract(t, closest);
	}

}
