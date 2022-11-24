package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractModule;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Lattice;
import fields.interfaces.Ring.FactorizationResult;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.FreeModule;
import fields.vectors.FreeSubModule;
import fields.vectors.Matrix;
import fields.vectors.Vector;
import util.Pair;

public class FractionalIdeal extends AbstractModule<IntE, NFE>
		implements Element<FractionalIdeal>, Lattice<NFE, Real, Vector<Real>> {
	private NumberField field;
	private NumberFieldIntegers order;
	private NumberFieldIdeal numerator;
	private NumberFieldIdeal denominator;
	private NumberFieldIdeal numeratorClearedDenominator;
	private IntE clearedDenominator;
	private List<NFE> basis;
	private Matrix<Fraction> toBasisBaseChange;
	private Matrix<Fraction> fromBasisBaseChange;
	private Matrix<Real> generatorsAsMatrix;

	public FractionalIdeal(NumberFieldIntegers order, NumberFieldIdeal numerator) {
		this(order, numerator, order.getUnitIdeal());
	}

	public FractionalIdeal(NumberFieldIntegers order, NumberFieldIdeal numerator, IntE denominator) {
		this(order, numerator, order.getIdeal(Collections.singletonList(order.getInteger(denominator))));
	}

	public FractionalIdeal(NumberFieldIntegers order, NumberFieldIdeal numerator, NumberFieldIdeal denominator) {
		this.order = order;
		this.field = order.numberField();
		if (numerator.equals(order.getZeroIdeal())) {
			this.numerator = order.getZeroIdeal();
			this.denominator = order.getUnitIdeal();
			return;
		}
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> numeratorFactors = order.idealFactorization(numerator);
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> denominatorFactors = order.idealFactorization(denominator);
		Set<Ideal<NFE>> primeIdeals = new TreeSet<>();
		primeIdeals.addAll(numeratorFactors.primeFactors());
		primeIdeals.addAll(denominatorFactors.primeFactors());
		Map<NumberFieldIdeal, Integer> numeratorReduced = new TreeMap<>();
		Map<NumberFieldIdeal, Integer> denominatorReduced = new TreeMap<>();
		for (Ideal<NFE> prime : primeIdeals) {
			int multiplicity = numeratorFactors.multiplicity(prime) - denominatorFactors.multiplicity(prime);
			if (multiplicity > 0) {
				numeratorReduced.put((NumberFieldIdeal) prime, multiplicity);
			}
			if (multiplicity < 0) {
				denominatorReduced.put((NumberFieldIdeal) prime, -multiplicity);
			}
		}
		this.numerator = order.fromFactorization(numeratorReduced);
		this.denominator = order.fromFactorization(denominatorReduced);
	}

	public static FractionalIdeal principalIdeal(NumberFieldIntegers order, NFE generator) {
		return new FractionalIdeal(order, order.getIdeal(Collections.singletonList(order.getNumerator(generator))),
				order.getIdeal(Collections.singletonList(order.getDenominator(generator))));
	}

	@Override
	public String toString() {
		return getNumerator() + "/" + getDenominator();
	}

	public NumberFieldIdeal getNumerator() {
		return numerator;
	}

	public NumberFieldIdeal getDenominator() {
		return denominator;
	}

	public boolean equals(Object o) {
		if (!(o instanceof FractionalIdeal)) {
			return false;
		}
		return compareTo((FractionalIdeal) o) == 0;
	}

	@Override
	public int compareTo(FractionalIdeal o) {
		int cmp = denominator.compareTo(o.denominator);
		if (cmp != 0) {
			return cmp;
		}
		return numerator.compareTo(o.numerator);
	}

	public Pair<NumberFieldIdeal, IntE> clearDenominator() {
		if (clearedDenominator == null) {
			Integers z = Integers.z();
			FactorizationResult<Ideal<NFE>, Ideal<NFE>> denominatorFactors = order.idealFactorization(denominator);
			Map<IntE, Integer> primes = new TreeMap<>();
			for (Ideal<NFE> primeFactor : denominatorFactors.primeFactors()) {
				NumberFieldIdeal ideal = (NumberFieldIdeal) primeFactor;
				if (!primes.containsKey(ideal.intGenerator())) {
					primes.put(ideal.intGenerator(), 0);
				}
				primes.put(ideal.intGenerator(),
						Math.max(primes.get(ideal.intGenerator()), denominatorFactors.multiplicity(primeFactor)));
			}
			IntE num = z.one();
			for (IntE prime : primes.keySet()) {
				num = z.multiply(num, z.power(prime, primes.get(prime)));
			}
			FractionalIdeal cleared = new FractionalIdeal(order,
					order.getIdeal(Collections.singletonList(order.getEmbedding(num))), denominator);
			if (!cleared.isInteger()) {
				throw new ArithmeticException("Could not clear the denominator");
			}
			numeratorClearedDenominator = order.multiply(this.numerator, cleared.numerator);
			clearedDenominator = num;
		}
		return new Pair<>(numeratorClearedDenominator, clearedDenominator);
	}

	public boolean isPrincipal() {
		return clearDenominator().getFirst().isPrincipal();
	}

	public boolean isInteger() {
		return denominator.equals(order.getUnitIdeal());
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public NFE zero() {
		return field.zero();
	}

	@Override
	public NFE add(NFE s1, NFE s2) {
		return field.add(s1, s2);
	}

	@Override
	public NFE negative(NFE s) {
		return field.negative(s);
	}

	@Override
	public NFE scalarMultiply(IntE t, NFE s) {
		return field.multiply(t, s);
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

	private List<Vector<IntE>> asVectors(List<NFE> s) {
		List<Vector<IntE>> result = new ArrayList<>();
		for (NFE t : s) {
			result.add(asVector(t));
		}
		return result;
	}

	@Override
	public boolean isLinearIndependent(List<NFE> s) {
		return new FreeModule<>(Integers.z(), field.degree()).isLinearIndependent(asVectors(s));
	}

	@Override
	public boolean isGeneratingModule(List<NFE> s) {
		return new FreeModule<>(Integers.z(), field.degree()).isGeneratingModule(asVectors(s));
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<NFE> s) {
		return new FreeModule<>(Integers.z(), field.degree()).nonTrivialCombinations(asVectors(s));
	}

	@Override
	public List<NFE> getModuleGenerators() {
		if (basis == null) {
			Pair<NumberFieldIdeal, IntE> cleared = clearDenominator();
			List<NFE> clearedGenerators = new ArrayList<>();
			for (NFE orderGenerator : order.getModuleGenerators()) {
				for (NFE numeratorGenerator : cleared.getFirst().generators()) {
					clearedGenerators.add(field.multiply(orderGenerator, numeratorGenerator));
				}
			}
			FreeSubModule<IntE, NFE> freeSubModule = new FreeSubModule<>(order, clearedGenerators);
			List<NFE> unreducedBasis = new ArrayList<>();
			for (NFE numeratorGenerator : freeSubModule.getModuleGenerators()) {
				unreducedBasis.add(field.divide(numeratorGenerator, field.getEmbedding(cleared.getSecond())));
			}
			this.basis = field.minkowskiEmbeddingSpace().latticeReduction(unreducedBasis, this);
			List<Vector<Fraction>> basisAsVectors = new ArrayList<>();
			for (NFE basisVector : basis) {
				basisAsVectors.add(field.asVector(basisVector));
			}
			this.fromBasisBaseChange = Matrix.fromColumns(basisAsVectors);
			this.toBasisBaseChange = field.matrixAlgebra().inverse(fromBasisBaseChange);
		}
		return basis;
	}

	@Override
	public Matrix<Real> generatorsAsMatrix() {
		if (generatorsAsMatrix == null) {
			List<Vector<Real>> asVectors = new ArrayList<>();
			for (NFE generator : getModuleGenerators()) {
				asVectors.add(embedding(generator));
			}
			generatorsAsMatrix = Matrix.fromColumns(asVectors);
		}
		return generatorsAsMatrix;
	}

	@Override
	public Vector<IntE> asVector(NFE s) {
		Vector<Fraction> inBasis = field.matrixAlgebra().multiply(toBasisBaseChange, field.asVector(s));
		List<IntE> result = new ArrayList<>();
		for (Fraction c : inBasis.asList()) {
			result.add(c.asInteger());
		}
		return new Vector<>(result);
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public NFE getRandomElement() {
		return fromVector(new FreeModule<IntE>(Integers.z(), field.degree()).getRandomElement());
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
	public Iterator<NFE> iterator() {
		throw new InfinityException();
	}

	@Override
	public FiniteRealVectorSpace getVectorSpace() {
		return field.minkowskiEmbeddingSpace();
	}

	@Override
	public Vector<Real> embedding(NFE t) {
		return field.minkowskiEmbedding(t);
	}

	@Override
	public int rank() {
		return field.degree();
	}
}
