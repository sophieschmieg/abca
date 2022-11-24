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
import fields.integers.FiniteRationalVectorSpace;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Integers.IntegerIdeal;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Lattice;
import fields.interfaces.MathMap;

public class EmbeddedRationalLattice extends AbstractModule<IntE, Vector<Fraction>>
		implements Lattice<Vector<Fraction>, Real, Vector<Real>> {
	private FiniteRealVectorSpace realSpace;
	private FiniteRationalVectorSpace rationalSpace;
	private List<Vector<Fraction>> basis;
	private Matrix<Real> realBaseChangeMatrix;
	private Matrix<Fraction> baseChangeMatrix;
	private MatrixModule<Fraction> matrixModule;
	private FreeModule<IntE> asIntModule;
	private MathMap<Fraction, Real> embeddingMap;

	public static EmbeddedRationalLattice fromRationalLattice(Reals r, RationalLattice lattice) {
		return new EmbeddedRationalLattice(r, lattice.getVectorSpace().dimension(), lattice.getModuleGenerators());
	}

	public static RationalLattice fromEmbeddedRationalLattice(EmbeddedRationalLattice lattice) {
		return new RationalLattice(new FiniteRationalVectorSpace(lattice.getVectorSpace().dimension()), lattice.getModuleGenerators());
	}

	public EmbeddedRationalLattice(Reals r, int dimension, List<Vector<Fraction>> generators) {
		this.rationalSpace = new FiniteRationalVectorSpace(dimension);
		this.realSpace = new FiniteRealVectorSpace(r, dimension);
		this.embeddingMap = new MathMap<>() {
			@Override
			public Real evaluate(Fraction t) {
				return r.getEmbedding(t);
			}};
		Matrix<Fraction> generatorMatrix = Matrix.fromColumns(generators);
		MatrixModule<Fraction> generatorMatrixModule = new MatrixModule<>(rationalSpace.getField(), dimension,
				generators.size());
		this.basis = new ArrayList<>();
		if (generatorMatrixModule.rank(generatorMatrix) < generators.size()) {
			Integers z = Integers.z();
			Rationals q = Rationals.q();
			IntE lcmDenominator = z.one();
			for (Vector<Fraction> generator : generators) {
				for (Fraction coeff : generator.asList()) {
					lcmDenominator = z.lcm(coeff.getDenominator(), lcmDenominator);
				}
			}
			IntE denominator = lcmDenominator;
			List<Vector<IntE>> integerList = new ArrayList<>();
			for (Vector<Fraction> generator : generators) {
				integerList.add(Vector.mapVector(new MathMap<Fraction, IntE>() {
					@Override
					public IntE evaluate(Fraction t) {
						return q.multiply(denominator, t).asInteger();
					}
				}, generator));
			}
			Matrix<IntE> asMatrix = Matrix.fromColumns(integerList);
			List<Vector<IntE>> reducedList = asMatrix.getModule(z).imageBasis(asMatrix);
			for (Vector<IntE> generator : reducedList) {
				basis.add(Vector.mapVector(new MathMap<IntE, Fraction>() {
					@Override
					public Fraction evaluate(IntE t) {
						return q.getFraction(t, denominator);
					}
				}, generator));
			}
		} else {
			this.basis.addAll(generators);
		}
		this.basis = realSpace.latticeReduction(this);
		this.baseChangeMatrix = Matrix.fromColumns(basis);
		this.realBaseChangeMatrix = Matrix.mapMatrix(embeddingMap, baseChangeMatrix);
		this.matrixModule = new MatrixModule<>(rationalSpace.getField(), rationalSpace.dimension(), basis.size());
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
	public Vector<Fraction> zero() {
		return rationalSpace.zero();
	}

	@Override
	public Vector<Fraction> add(Vector<Fraction> s1, Vector<Fraction> s2) {
		return rationalSpace.add(s1, s2);
	}

	@Override
	public Vector<Fraction> negative(Vector<Fraction> s) {
		return rationalSpace.negative(s);
	}

	@Override
	public Vector<Fraction> scalarMultiply(IntE t, Vector<Fraction> s) {
		return rationalSpace.scalarMultiply(t.getValue(), s);
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
	public List<Vector<IntE>> getSyzygies() {
		return Collections.emptyList();
	}

	@Override
	public boolean isLinearIndependent(List<Vector<Fraction>> s) {
		return rationalSpace.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<Vector<Fraction>> s) {
		return asIntModule.isGeneratingModule(asIntVectors(s));
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<Vector<Fraction>> s) {
		return asIntModule.nonTrivialCombinations(asIntVectors(s));
	}

	private List<Vector<IntE>> asIntVectors(List<Vector<Fraction>> s) {
		List<Vector<IntE>> result = new ArrayList<>();
		for (Vector<Fraction> vector : s) {
			result.add(asVector(vector));
		}
		return result;
	}

	@Override
	public List<Vector<Fraction>> getModuleGenerators() {
		return basis;
	}

	@Override
	public Matrix<Real> generatorsAsMatrix() {
		return realBaseChangeMatrix;
	}

	@Override
	public Vector<IntE> asVector(Vector<Fraction> s) {
		Vector<Fraction> realSolution = matrixModule.solve(baseChangeMatrix, s);
		List<IntE> integerSolution = new ArrayList<>();
		for (Fraction realCoefficient : realSolution.asList()) {
			integerSolution.add(realCoefficient.round());
		}
		return new Vector<>(integerSolution);
	}

	public boolean contains(Vector<Fraction> vector) {
		if (!matrixModule.hasSolution(baseChangeMatrix, vector)) {
			return false;
		}
		Vector<Fraction> solution = matrixModule.solve(baseChangeMatrix, vector);
		for (Fraction coefficient : solution.asList()) {
			if (!coefficient.getDenominator().equals(Integers.z().one())) {
				return false;
			}
		}
		return true;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public Vector<Fraction> getRandomElement() {
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
	public Iterator<Vector<Fraction>> iterator() {
		throw new InfinityException();
	}

	@Override
	public FiniteRealVectorSpace getVectorSpace() {
		return realSpace;
	}

	@Override
	public Vector<Real> embedding(Vector<Fraction> t) {
		return Vector.mapVector(embeddingMap, t);
	}

}
