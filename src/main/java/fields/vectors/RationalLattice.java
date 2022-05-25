package fields.vectors;

import java.math.BigInteger;
import java.util.ArrayList;
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
import fields.interfaces.InnerProductSpace;
import fields.interfaces.Lattice;
import fields.interfaces.RealInnerProductSpace;

public class RationalLattice extends AbstractModule<IntE, Vector<Real>>
		implements Lattice<Vector<Real>, Real, Vector<Real>> {
	private RealInnerProductSpace<Real, Vector<Real>> space;
	private List<Vector<Real>> basis;
	private Matrix<Real> baseChangeMatrix;
	private MatrixModule<Real> matrixModule;
	private FreeModule<IntE> asIntModule;

	public RationalLattice(FiniteRealVectorSpace space, List<Vector<Real>> generators) {
		this(space, generators, 0.75);
	}

	public RationalLattice(FiniteRealVectorSpace space, List<Vector<Real>> generators, double delta) {
		this.space = space;
		Matrix<Real> generatorMatrix = Matrix.fromColumns(generators);
		MatrixModule<Real> generatorMatrixModule = new MatrixModule<>(space.getField(), space.dimension(),
				generators.size());
		this.basis = new ArrayList<>();
		if (generatorMatrixModule.rank(generatorMatrix) < generators.size()) {
			Integers z = Integers.z();
			Rationals q = Rationals.q();
			Reals r = space.getValueField().getReals();
			List<Vector<Real>> independentSubset = space.linearIndependentSubSet(generators);
			int[] subSetIndeces = new int[independentSubset.size()];
			for (int i = 0; i < subSetIndeces.length; i++) {
				for (int j = 0; j < generators.size(); j++) {
					Vector<Real> generator = generators.get(j);
					if (independentSubset.get(i).equals(generator)) {
						subSetIndeces[i] = j;
						break;
					}
				}
			}
			Matrix<Real> baseChange = Matrix.fromColumns(independentSubset);
			List<List<IntE>> integerKernel = new ArrayList<>();
			for (int i = 0; i < generators.size(); i++) {
				Vector<Real> generator = generators.get(i);
				Vector<Real> inBasis = generatorMatrixModule.solve(baseChange, generator);
				List<Fraction> asFractions = new ArrayList<>();
				IntE lcm = z.one();
				for (Real coefficient : inBasis.asList()) {
					Fraction asFraction = r.roundToFraction(coefficient, r.precision() / 2);
					asFractions.add(asFraction);
					lcm = z.lcm(asFraction.getDenominator(), lcm);
				}
				List<IntE> integerRelation = new ArrayList<>();
				for (int j = 0; j < generators.size(); j++) {
					integerRelation.add(z.zero());
				}
				integerRelation.set(i, z.negative(lcm));
				for (int j = 0; j < independentSubset.size(); j++) {
					integerRelation.set(subSetIndeces[j], z.add(integerRelation.get(subSetIndeces[j]),
							q.multiply(lcm, asFractions.get(j)).asInteger()));
				}
				integerKernel.add(integerRelation);
			}
//			List<Vector<Real>> kernel = generatorMatrixModule.kernelBasis(generatorMatrix);
//			for (Vector<Real> kernelVector : kernel) {
//				Real coeff = null;
//				for (int i = 0; i < generators.size(); i++) {
//					coeff = kernelVector.get(i + 1);
//					if (r.abs(coeff).compareTo(r.getPowerOfTwo(-100)) > 0) {
//						break;
//					}
//				}
//				kernelVector = generatorMatrixModule.domain().scalarMultiply(r.inverse(coeff), kernelVector);
//				IntE lcm = z.one();
//				List<Fraction> asFractions = new ArrayList<>();
//				for (Real coefficient : kernelVector.asList()) {
//					Fraction asFraction = r.roundToFraction(coefficient, 100);
//					asFractions.add(asFraction);
//					lcm = z.lcm(asFraction.getDenominator(), lcm);
//				}
//				List<IntE> integerVector = new ArrayList<>();
//				for (Fraction coefficient : asFractions) {
//					integerVector.add(q.multiply(lcm, coefficient).asInteger());
//				}
//				integerKernel.add(integerVector);
//			}
			FreeModule<IntE> freeModule = new FreeModule<>(z, generators.size());
			GenericPIDModule<IntE, Vector<IntE>> modKernel = GenericPIDModule.fromSyzygies(freeModule, integerKernel);
			List<IntE> diagonalRanks = modKernel.diagonalRanks();
			for (int i = 0; i < generators.size(); i++) {
				IntE diagonalRank = diagonalRanks.get(i);
				if (!diagonalRank.equals(z.zero())) {
					continue;
				}
				Vector<IntE> diagonal = modKernel.lift(modKernel.fromDiagonalVector(freeModule.getUnitVector(i + 1)));
				Vector<Real> basisVector = space.zero();
				for (int j = 0; j < generators.size(); j++) {
					basisVector = space.add(space.scalarMultiply(diagonal.get(j + 1).getValue(), generators.get(j)),
							basisVector);
				}
				this.basis.add(basisVector);
			}
		} else {
			this.basis.addAll(generators);
		}
		this.basis = space.latticeReduction(this, delta);
		this.baseChangeMatrix = Matrix.fromColumns(basis);
		this.matrixModule = new MatrixModule<>(space.getField(), space.dimension(), basis.size());
		this.asIntModule = new FreeModule<>(Integers.z(), basis.size());
	}

	@Override
	public String toString() {
		return basis.toString();
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public Vector<Real> zero() {
		return space.zero();
	}

	@Override
	public Vector<Real> add(Vector<Real> s1, Vector<Real> s2) {
		return space.add(s1, s2);
	}

	@Override
	public Vector<Real> negative(Vector<Real> s) {
		return space.negative(s);
	}

	@Override
	public Vector<Real> scalarMultiply(IntE t, Vector<Real> s) {
		return space.scalarMultiply(t.getValue(), s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public int rank() {
		return basis.size();
	}

	@Override
	public IntegerIdeal annihilator() {
		return Integers.z().getZeroIdeal();
	}

	@Override
	public boolean isLinearIndependent(List<Vector<Real>> s) {
		return space.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<Vector<Real>> s) {
		return asIntModule.isGeneratingModule(asIntVectors(s));
	}

	@Override
	public List<List<IntE>> nonTrivialCombinations(List<Vector<Real>> s) {
		return asIntModule.nonTrivialCombinations(asIntVectors(s));
	}

	private List<Vector<IntE>> asIntVectors(List<Vector<Real>> s) {
		List<Vector<IntE>> result = new ArrayList<>();
		for (Vector<Real> vector : s) {
			result.add(asVector(vector));
		}
		return result;
	}

	@Override
	public List<Vector<Real>> getModuleGenerators() {
		return basis;
	}

	@Override
	public Matrix<Real> generatorsAsMatrix() {
		return baseChangeMatrix;
	}

	@Override
	public Vector<IntE> asVector(Vector<Real> s) {
		Vector<Real> realSolution = matrixModule.solve(baseChangeMatrix, s);
		List<IntE> integerSolution = new ArrayList<>();
		for (Real realCoefficient : realSolution.asList()) {
			integerSolution.add(realCoefficient.round());
		}
		return new Vector<>(integerSolution);
	}

	public boolean contains(Vector<Real> vector) {
		if (!matrixModule.hasSolution(baseChangeMatrix, vector)) {
			return false;
		}
		Vector<Real> realSolution = matrixModule.solve(baseChangeMatrix, vector);
		Reals r = (Reals) space.getField();
		Integers z = Integers.z();
		for (Real coefficient : realSolution.asList()) {
			if (!r.roundToFraction(coefficient, 100).getDenominator().equals(z.one())) {
				return false;
			}
		}
		return true;
	}

	@Override
	public Exactness exactness() {
		return Exactness.FLOATING_POINT;
	}

	@Override
	public Vector<Real> getRandomElement() {
		return fromVector(asIntModule.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<Vector<Real>> iterator() {
		throw new InfinityException();
	}

	@Override
	public InnerProductSpace<Real, Vector<Real>> getVectorSpace() {
		return space;
	}

	@Override
	public Vector<Real> embedding(Vector<Real> t) {
		return t;
	}

}
