package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractIdeal;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Algebra;
import fields.interfaces.DedekindRing;
import fields.interfaces.Ideal;
import fields.interfaces.LocalRing;
import fields.interfaces.LocalRing.OkutsuType;
import fields.interfaces.LocalRing.TheMontesResult;
import fields.interfaces.MathMap;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.Value;
import fields.numberfields.NumberField.NFE;
import fields.vectors.FreeModule;
import fields.vectors.FreeSubModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.Identity;
import util.MiscAlgorithms;
import util.SingletonSortedMap;

public class NumberFieldIntegers extends AbstractAlgebra<IntE, NFE>
		implements Algebra<IntE, NFE>, DedekindRing<NFE, NFE, FFE> {
	private NumberField field;
	private UnivariatePolynomial<IntE> minimalPolynomial;
	private List<NFE> integralBasis;
	private Matrix<Fraction> toIntegralBasis;
	private Matrix<Fraction> fromIntegralBasis;
	private Map<IntE, List<NumberFieldIdeal>> idealsOverPrime;
	private Map<FactorizationResult<Ideal<NFE>, Ideal<NFE>>, NumberFieldIdeal> idealsByFactorization;
	private NumberFieldIdeal unitIdeal;
	private NumberFieldIdeal zeroIdeal;
	private FreeModule<IntE> asFreeModule;
	private Map<Ideal<NFE>, LocalizedNumberField> localizations;

	NumberFieldIntegers(NumberField field) {
		this.field = field;
		this.idealsOverPrime = new TreeMap<>();
		this.idealsByFactorization = new TreeMap<>();
		this.asFreeModule = new FreeModule<>(Integers.z(), field.degree());
		this.localizations = new TreeMap<>();
	}

	@Override
	public String toString() {
		return "O_" + field.toString();
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public NFE getRandomElement() {
		return fromVector(asFreeModule.getRandomElement());
	}

	@Override
	public boolean isCommutative() {
		return true;
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
	public NFE zero() {
		return field.zero();
	}

	@Override
	public NFE one() {
		return field.one();
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public NFE add(NFE t1, NFE t2) {
		return field.add(t1, t2);
	}

	@Override
	public NFE negative(NFE t) {
		return field.negative(t);
	}

	@Override
	public NFE multiply(NFE t1, NFE t2) {
		return field.multiply(t1, t2);
	}

	@Override
	public boolean isUnit(NFE t) {
		return Integers.z().isUnit(field.norm(t).getNumerator());
	}

	@Override
	public NFE inverse(NFE t) {
		if (!isUnit(t)) {
			throw new ArithmeticException(t + " is not invertible");
		}
		return field.inverse(t);
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isZeroDivisor(NFE t) {
		return false;
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	@Override
	public FactorizationResult<NFE, NFE> uniqueFactorization(NFE t) {
		throw new ArithmeticException("Not a UFD!");
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		return true;
	}

	@Override
	public boolean isDivisible(NFE dividend, NFE divisor) {
		return isElement(field.divide(dividend, divisor));
	}

	@Override
	public QuotientAndRemainderResult<NFE> quotientAndRemainder(NFE dividend, NFE divisor) {
		Ideal<NFE> divisorIdeal = getIdeal(Collections.singletonList(divisor));
		return new QuotientAndRemainderResult<>(divisorIdeal.generate(dividend).get(0), divisorIdeal.residue(dividend));
	}

	@Override
	public BigInteger euclidMeasure(NFE t) {
		throw new ArithmeticException("Not a Euclidean ring!");
	}

	@Override
	public NFE projectToUnit(NFE t) {
		return one();
	}

	@Override
	public Iterable<NFE> getUnits() {
		if (field.degree() > 2) {
			throw new UnsupportedOperationException("Not implemented!");
		}
		if (field.degree() == 1) {
			List<NFE> units = new ArrayList<>();
			units.add(getInteger(-1));
			units.add(getInteger(1));
			return units;
		}
		Integers z = Integers.z();
		if (field.discriminant().compareTo(z.zero()) < 0) {
			List<NFE> units = new ArrayList<>();
			units.add(getInteger(-1));
			units.add(getInteger(1));
			if (field.discriminant().equals(z.negative(z.getInteger(4)))) {
				units.addAll(field.sqrt(field.getInteger(-1)).keySet());
			} else if (field.discriminant().equals(z.negative(z.getInteger(3)))) {
				UnivariatePolynomialRing<NFE> polynomials = field.getUnivariatePolynomialRing();
				UnivariatePolynomial<NFE> eisenstein = polynomials.getPolynomial(field.one(), field.one(), field.one());
				units.addAll(field.roots(eisenstein).keySet());
			}
			return units;
		}
		IntE discriminant = field.discriminant();
		boolean mod23;
		if (discriminant.getValue().mod(BigInteger.valueOf(4)).equals(BigInteger.ZERO)) {
			discriminant = Integers.z().divideChecked(discriminant, z.getInteger(4));
			mod23 = true;
		} else {
			mod23 = false;
		}
		List<Vector<IntE>> pells = MiscAlgorithms.pellsEquation(discriminant, z.getInteger(mod23 ? -1 : -4), false);
		pells.addAll(MiscAlgorithms.pellsEquation(discriminant, z.getInteger(mod23 ? 1 : 4), false));
		Vector<IntE> pell = pells.get(0);
		NFE alpha = field.sqrt(field.getInteger(discriminant)).keySet().iterator().next();
		final NFE fundamental;
		if (mod23) {
			fundamental = add(getInteger(pell.get(1)), multiply(pell.get(2), alpha));
		} else {
			fundamental = divide(add(getInteger(pell.get(1)), multiply(pell.get(2), alpha)), getInteger(2));
		}
		return new Iterable<NumberField.NFE>() {

			@Override
			public Iterator<NFE> iterator() {
				return new Iterator<>() {
					Iterator<IntE> it = z.iterator();
					IntE current = null;

					@Override
					public boolean hasNext() {
						return true;
					}

					@Override
					public NFE next() {
						if (current == null) {
							current = it.next();
							return power(fundamental, current.intValueExact());
						}
						NFE result = negative(power(fundamental, current.intValueExact()));
						current = null;
						return result;
					}
				};
			}
		};
	}

	@Override
	public int krullDimension() {
		return 0;
	}

	@Override
	public IdealResult<NFE, NumberFieldIdeal> getIdealWithTransforms(List<NFE> generators) {
		NumberFieldIdeal ideal = getIdeal(generators);
		List<NFE> integralBasis = getModuleGenerators();
		List<Vector<IntE>> columns = new ArrayList<>();
		for (NFE generator : generators) {
			for (NFE basisElement : integralBasis) {
				NFE multiplied = multiply(basisElement, generator);
				columns.add(asVector(multiplied));
			}
		}
		Matrix<IntE> matrix = Matrix.fromColumns(columns);
		MatrixModule<IntE> mm = new MatrixModule<>(Integers.z(), integralBasis.size(), columns.size());
		List<List<NFE>> transforms = new ArrayList<>();
		for (NFE basisElement : ideal.generators()) {
			Vector<IntE> rhs = asVector(basisElement);
			List<IntE> solution = mm.solve(matrix, rhs).asList();
			List<NFE> transform = new ArrayList<>();
			for (int i = 0; i < generators.size(); i++) {
				transform.add(fromVector(
						new Vector<>(solution.subList(i * integralBasis.size(), (i + 1) * integralBasis.size()))));
			}
			transforms.add(transform);
		}
		return new IdealResult<>(transforms, generators, ideal);
	}

	@Override
	public NumberFieldIdeal getIdeal(List<NFE> generators) {
		Integers z = Integers.z();
		IntE normGcd = z.zero();
		for (NFE generator : generators) {
			normGcd = z.gcd(normGcd, field.norm(generator).asInteger());
		}
		if (normGcd.equals(z.zero())) {
			List<NFE> expression = new ArrayList<>();
			for (int i = 0; i < generators.size(); i++) {
				expression.add(zero());
			}
			return getZeroIdeal();
		}
		FactorizationResult<IntE, IntE> factorization = z.uniqueFactorization(normGcd);
		return getIdealWithFactorization(generators, factorization);
	}

	public NumberFieldIdeal getIdealIfSmoothOver(List<NFE> generators, Set<IntE> smoothnessBase) {
		Integers z = Integers.z();
		IntE normGcd = z.zero();
		for (NFE generator : generators) {
			normGcd = z.gcd(normGcd, field.norm(generator).asInteger());
		}
		if (normGcd.equals(z.zero())) {
			return getZeroIdeal();
		}
		SortedMap<IntE, Integer> factors = new TreeMap<>();
		for (IntE prime : smoothnessBase) {
			int power = 0;
			while (z.isDivisible(normGcd, prime)) {
				normGcd = z.divide(normGcd, prime);
				power++;
			}
			if (power != 0) {
				factors.put(prime, power);
			}
		}
		if (!z.isUnit(normGcd)) {
			return null;
		}
		return getIdealWithFactorization(generators, new FactorizationResult<>(normGcd, factors));
	}

	private NumberFieldIdeal getIdealWithFactorization(List<NFE> generators,
			FactorizationResult<IntE, IntE> normFactorization) {
		if (normFactorization.primeFactors().isEmpty()) {
			return getUnitIdeal();
		}
		Integers z = Integers.z();
		SortedMap<IntE, FactorizationResult<Ideal<NFE>, Ideal<NFE>>> factorizations = new TreeMap<>();
		SortedMap<Ideal<NFE>, Integer> allIdealFactors = new TreeMap<>();
		for (IntE prime : normFactorization.primeFactors()) {
			SortedMap<Ideal<NFE>, Integer> idealFactors = new TreeMap<>();
			List<NumberFieldIdeal> ideals = idealsOver(z.getIdeal(Collections.singletonList(prime)));
			for (NumberFieldIdeal ideal : ideals) {
				Value value = Value.INFINITY;
				for (NFE generator : generators) {
					if (generator.equals(zero())) {
						continue;
					}
					value = value.min(ideal.maximumPowerContains(generator));
				}
				if (!value.equals(Value.ZERO)) {
					idealFactors.put(ideal, value.value());
					allIdealFactors.put(ideal, value.value());
				}
			}
			if (normFactorization.primeFactors().size() == 1 && idealFactors.size() == 1
					&& idealFactors.get(idealFactors.firstKey()) == 1) {
				return (NumberFieldIdeal) idealFactors.firstKey();
			}
			factorizations.put(prime, new FactorizationResult<>(getUnitIdeal(), idealFactors));
		}
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = new FactorizationResult<>(getUnitIdeal(), allIdealFactors);
		if (!idealsByFactorization.containsKey(factorization)) {
			idealsByFactorization.put(factorization, new NumberFieldIdeal(factorizations));
		}
		return idealsByFactorization.get(factorization);
	}

	@Override
	public NumberFieldIdeal multiply(Ideal<NFE> t1, Ideal<NFE> t2) {
		NumberFieldIdeal ideal1 = (NumberFieldIdeal) t1;
		NumberFieldIdeal ideal2 = (NumberFieldIdeal) t2;
		if (ideal1.zeroIdeal || ideal2.zeroIdeal) {
			return getZeroIdeal();
		}
		Set<Ideal<NFE>> primeFactors = new TreeSet<>();
		primeFactors.addAll(ideal1.idealFactorization.primeFactors());
		primeFactors.addAll(ideal2.idealFactorization.primeFactors());
		Map<NumberFieldIdeal, Integer> result = new TreeMap<>();
		for (Ideal<NFE> ideal : primeFactors) {
			result.put((NumberFieldIdeal) ideal,
					ideal1.idealFactorization.multiplicity(ideal) + ideal2.idealFactorization.multiplicity(ideal));
		}
		return fromFactorization(result);
	}

	@Override
	public Ideal<NFE> intersect(Ideal<NFE> t1, Ideal<NFE> t2) {
		NumberFieldIdeal ideal1 = (NumberFieldIdeal) t1;
		NumberFieldIdeal ideal2 = (NumberFieldIdeal) t2;
		if (ideal1.zeroIdeal || ideal2.zeroIdeal) {
			return getZeroIdeal();
		}
		Set<Ideal<NFE>> primeFactors = new TreeSet<>();
		primeFactors.addAll(ideal1.idealFactorization.primeFactors());
		primeFactors.addAll(ideal2.idealFactorization.primeFactors());
		Map<NumberFieldIdeal, Integer> result = new TreeMap<>();
		for (Ideal<NFE> ideal : primeFactors) {
			result.put((NumberFieldIdeal) ideal, Math.max(ideal1.idealFactorization.multiplicity(ideal),
					ideal2.idealFactorization.multiplicity(ideal)));
		}
		return fromFactorization(result);
	}

	@Override
	public NumberFieldIdeal radical(Ideal<NFE> t) {
		NumberFieldIdeal ideal = (NumberFieldIdeal) t;
		if (ideal.zeroIdeal) {
			return getZeroIdeal();
		}
		Map<NumberFieldIdeal, Integer> result = new TreeMap<>();
		for (Ideal<NFE> idealFactor : ideal.idealFactorization.primeFactors()) {
			result.put((NumberFieldIdeal) idealFactor, 1);
		}
		return fromFactorization(result);
	}

	@Override
	public NumberFieldIdeal power(Ideal<NFE> t, int power) {
		Map<NumberFieldIdeal, Integer> factors = new TreeMap<>();
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> originalFactors = idealFactorization(t);
		for (Ideal<NFE> prime : originalFactors.primeFactors()) {
			factors.put((NumberFieldIdeal) prime, power * originalFactors.multiplicity(prime));
		}
		return fromFactorization(factors);
	}

	NumberFieldIdeal fromFactorization(Map<NumberFieldIdeal, Integer> factors) {
		if (factors.size() == 0) {
			return getUnitIdeal();
		}
		if (factors.size() == 1 && factors.get(factors.keySet().iterator().next()) == 1) {
			return factors.keySet().iterator().next();
		}
		SortedMap<Ideal<NFE>, Integer> cast = new TreeMap<>();
		Map<IntE, SortedMap<Ideal<NFE>, Integer>> perPrimeFactorization = new TreeMap<>();
		for (Ideal<NFE> ideal : factors.keySet()) {
			NumberFieldIdeal primeIdeal = (NumberFieldIdeal) ideal;
			IntE prime = primeIdeal.intGenerator;
			if (!perPrimeFactorization.containsKey(prime)) {
				perPrimeFactorization.put(prime, new TreeMap<>());
			}
			perPrimeFactorization.get(prime).put(ideal, factors.get(ideal));
			cast.put(ideal, factors.get(ideal));
		}
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> asFactorizationResult = new FactorizationResult<>(getUnitIdeal(), cast);
		if (!idealsByFactorization.containsKey(asFactorizationResult)) {
			Map<IntE, FactorizationResult<Ideal<NFE>, Ideal<NFE>>> perPrime = new TreeMap<>();
			for (IntE prime : perPrimeFactorization.keySet()) {
				perPrime.put(prime, new FactorizationResult<>(getUnitIdeal(), perPrimeFactorization.get(prime)));
			}
			idealsByFactorization.put(asFactorizationResult, new NumberFieldIdeal(perPrime));
		}
		return idealsByFactorization.get(asFactorizationResult);
	}

	@Override
	public FactorizationResult<Ideal<NFE>, Ideal<NFE>> idealFactorization(Ideal<NFE> t) {
		NumberFieldIdeal ideal = (NumberFieldIdeal) t;
		return ideal.idealFactorization;
	}

	@Override
	public NumberFieldIdeal getUnitIdeal() {
		if (unitIdeal == null) {
			this.unitIdeal = new NumberFieldIdeal(false);
		}
		return unitIdeal;
	}

	@Override
	public NumberFieldIdeal getZeroIdeal() {
		if (zeroIdeal == null) {
			this.zeroIdeal = new NumberFieldIdeal(true);
		}
		return zeroIdeal;
	}

	public List<NumberFieldIdeal> idealsOver(Ideal<IntE> integerIdeal) {
		IntE generator = integerIdeal.generators().get(0);
		FactorizationResult<IntE, IntE> primes = Integers.z().uniqueFactorization(generator);
		if (!primes.isIrreducible()) {
			throw new ArithmeticException("Expected prime ideal");
		}
		IntE prime = primes.firstPrimeFactor();
		if (idealsOverPrime.containsKey(prime)) {
			return idealsOverPrime.get(prime);
		}
		LocalRing<Fraction, PFE> zp = Integers.z().localize(prime);
		PrimeField fp = PrimeField.getPrimeField(prime);
		Rationals q = Rationals.q();
		UnivariatePolynomial<Fraction> zpMinimalPolynomial = q.getUnivariatePolynomialRing()
				.getEmbedding(minimalPolynomial(), new MathMap<>() {
					@Override
					public Fraction evaluate(IntE t) {
						return q.getInteger(t);
					}
				});
		TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes = zp.theMontesAlgorithm(zpMinimalPolynomial,
				fp.getExtension(fp.getUnivariatePolynomialRing().getVar()));
		List<NumberFieldIdeal> ideals = new ArrayList<>();
		List<List<UnivariatePolynomial<Fraction>>> integralBasis = zp.integralBasis(zpMinimalPolynomial, theMontes,
				true);
		for (int i = 0; i < theMontes.getTypes().size(); i++) {
			ideals.add(new NumberFieldIdeal(zp, theMontes.getTypes(), i, integralBasis));
		}
		idealsOverPrime.put(prime, ideals);
		return ideals;
	}

	private UnivariatePolynomial<IntE> minimalPolynomial() {
		if (this.minimalPolynomial == null) {
			UnivariatePolynomial<Fraction> minimalPolynomial = field.minimalPolynomial();
			Integers z = Integers.z();
			Rationals q = Rationals.q();
			UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
			if (!minimalPolynomial.leadingCoefficient().equals(q.one())) {
				throw new ArithmeticException("Minimal polynomial should be normalized");
			}
			IntE denom = z.one();
			for (int i = 0; i < minimalPolynomial.degree(); i++) {
				denom = z.lcm(denom, minimalPolynomial.univariateCoefficient(i).getDenominator());
			}
			minimalPolynomial = rationalPolynomials
					.normalize(rationalPolynomials.substitute(minimalPolynomial, Collections.singletonList(
							rationalPolynomials.divideScalar(rationalPolynomials.getVar(), q.getInteger(denom)))));
			this.minimalPolynomial = z.getUnivariatePolynomialRing().getEmbedding(minimalPolynomial, new MathMap<>() {
				@Override
				public IntE evaluate(Fraction t) {
					return t.getNumerator();
				}
			});
		}
		return minimalPolynomial;
	}

	IntE discriminant(List<NFE> generators) {
		if (!field.isBasis(generators)) {
			throw new ArithmeticException("Discriminant needs a basis");
		}
		List<List<Fraction>> traceForm = new ArrayList<>();
		for (int i = 0; i < generators.size(); i++) {
			if (!field.isInteger(generators.get(i))) {
				throw new ArithmeticException("Not an integer element");
			}
			List<Fraction> traceFormRow = new ArrayList<>();
			for (int j = 0; j < generators.size(); j++) {
				traceFormRow.add(field.trace(field.multiply(generators.get(i), generators.get(j))));
			}
			traceForm.add(traceFormRow);
		}
		Matrix<Fraction> trace = new Matrix<>(traceForm);
		return field.matrixAlgebra().determinant(trace).asInteger();
	}

	private Set<IntE> potentiallyRamifiedPrimes() {
		Integers z = Integers.z();
		return z.uniqueFactorization(z.getUnivariatePolynomialRing().discriminant(minimalPolynomial())).primeFactors();
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
	public boolean isLinearIndependent(List<NFE> s) {
		List<Vector<IntE>> asVectors = new ArrayList<>();
		for (NFE t : s) {
			asVectors.add(asVector(t));
		}
		return asFreeModule.isLinearIndependent(asVectors);
	}

	@Override
	public boolean isGeneratingModule(List<NFE> s) {
		List<Vector<IntE>> asVectors = new ArrayList<>();
		for (NFE t : s) {
			asVectors.add(asVector(t));
		}
		return asFreeModule.isGeneratingModule(asVectors);
	}

	@Override
	public List<List<IntE>> nonTrivialCombinations(List<NFE> s) {
		List<Vector<IntE>> asVectors = new ArrayList<>();
		for (NFE t : s) {
			asVectors.add(asVector(t));
		}
		return asFreeModule.nonTrivialCombinations(asVectors);
	}

	public boolean isElement(NFE t) {
		Rationals q = Rationals.q();
		if (toIntegralBasis == null) {
			UnivariatePolynomial<Fraction> minPoly = field.minimalPolynomial(t);
			if (!minPoly.leadingCoefficient().equals(q.one())) {
				minPoly = q.getUnivariatePolynomialRing().normalize(minPoly);
			}
			for (int i = 0; i < minPoly.degree(); i++) {
				if (!q.isInteger(minPoly.univariateCoefficient(i))) {
					return false;
				}
			}
			return true;
		}
		Vector<Fraction> integerVector = field.matrixAlgebra().multiply(toIntegralBasis, field.asVector(t));
		for (Fraction c : integerVector.asList()) {
			if (!q.isInteger(c)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public List<NFE> getModuleGenerators() {
		if (integralBasis == null) {
			Rationals q = Rationals.q();
			Integers z = Integers.z();
			if (field.degree() == 1) {
				this.integralBasis = Collections.singletonList(one());
				this.fromIntegralBasis = field.matrixAlgebra().one();
				this.toIntegralBasis = field.matrixAlgebra().one();
				return integralBasis;
			}
			UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
			List<IntE> primes = new ArrayList<>();
			primes.addAll(potentiallyRamifiedPrimes());
			Map<IntE, List<UnivariatePolynomial<Fraction>>> integralBasisPerPrime = new TreeMap<>();
			Map<IntE, List<Integer>> valuesPerPrime = new TreeMap<>();
			BigInteger allPrimes = BigInteger.ONE;
			UnivariatePolynomial<Fraction> zpMinimalPolynomial = q.getUnivariatePolynomialRing()
					.getEmbedding(minimalPolynomial(), new MathMap<>() {
						@Override
						public Fraction evaluate(IntE t) {
							return q.getInteger(t);
						}
					});
			for (IntE prime : primes) {
				allPrimes = allPrimes.multiply(prime.getValue());
				PrimeField fp = PrimeField.getPrimeField(prime.getValue());
				LocalRing<Fraction, PFE> zp = z.localize(prime.getValue());
				TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes = zp.theMontesAlgorithm(
						zpMinimalPolynomial, fp.getExtension(fp.getUnivariatePolynomialRing().getVar()));
				List<UnivariatePolynomial<Fraction>> integral = zp.triagonalizeIntegralBasis(zpMinimalPolynomial,
						zp.integralBasis(zpMinimalPolynomial, theMontes, false));
				List<Integer> valueList = new ArrayList<>();
				for (int i = 0; i < minimalPolynomial().degree(); i++) {
					valueList.add(-zp.fieldOfFractions().valuation(integral.get(i).leadingCoefficient()).value());
				}
				integralBasisPerPrime.put(prime, integral);
				valuesPerPrime.put(prime, valueList);
			}
			this.integralBasis = new ArrayList<>();
			for (int i = 0; i < minimalPolynomial().degree(); i++) {
				UnivariatePolynomial<Fraction> b = polynomials.zero();
				IntE primeProduct = z.one();
				IntE primeProductDivisor = z.one();
				for (IntE prime : primes) {
					IntE primePower = z.power(prime, valuesPerPrime.get(prime).get(i));
					primeProductDivisor = z.multiply(primePower, primeProductDivisor);
					primeProduct = z.multiply(prime, primePower, primeProduct);
				}
				for (IntE prime : primes) {
					BigInteger primePower = z.power(prime, valuesPerPrime.get(prime).get(i) + 1).getValue();
					IntE allOthers = z.divideChecked(primeProduct, z.getInteger(prime));
					IntE multiplier = z.multiply(primeProduct.getValue().divide(primePower).modInverse(primePower),
							allOthers);
					b = polynomials.add(b, polynomials.multiply(multiplier, integralBasisPerPrime.get(prime).get(i)));
				}
				BigInteger product = primeProduct.getValue();
				b = polynomials.getEmbedding(b, new MathMap<>() {
					@Override
					public Fraction evaluate(Fraction t) {
						BigInteger rounded = t.getNumerator().getValue().divide(t.getDenominator().getValue());
						BigInteger quotient = rounded.divide(product);
						return q.subtract(t, q.getInteger(quotient.multiply(product)));
					}
				});
				b = polynomials.divideScalar(b, q.getEmbedding(primeProductDivisor));
				if (b.degree() != i) {
					throw new ArithmeticException("expected integral basis to be triagonal!");
				}
//				for (int j = b.degree() - 1; j >= 0; j--) {
//					Fraction c = b.univariateCoefficient(j);
//					Fraction quotient = q.divide(c, polynomialBasis.get(j).leadingCoefficient());
//					IntE intQuotient = z.getInteger(
//							quotient.getNumerator().getValue().divide(quotient.getDenominator().getValue()));
//					b = polynomials.toUnivariate(
//							polynomials.subtract(b, polynomials.multiply(intQuotient, polynomialBasis.get(j))));
//				}
				NFE element = field.fromPolynomial(b);
				integralBasis.add(element);
//				asVectors.add(field.asVector(element));
//				polynomialBasis.add(b);
			}
			this.integralBasis = field.latticeReduction(integralBasis);
			List<Vector<Fraction>> asVectors = new ArrayList<>();
			for (NFE b : integralBasis) {
				asVectors.add(field.asVector(b));
			}
			this.fromIntegralBasis = Matrix.fromColumns(asVectors);
			this.toIntegralBasis = field.matrixAlgebra().inverse(fromIntegralBasis);
		}
		return integralBasis;
	}

	@Override
	public Vector<IntE> asVector(NFE s) {
		Vector<Fraction> asFieldVector = field.asVector(s);
		Vector<Fraction> asIntVector = field.matrixAlgebra().multiply(toIntegralBasisBaseChange(), asFieldVector);
		List<IntE> ints = new ArrayList<>();
		for (Fraction c : asIntVector.asList()) {
			ints.add(c.asInteger());
		}
		return new Vector<>(ints);
	}

	public NFE roundToInteger(NFE t) {
		Vector<Fraction> inIntegralBasis = field.matrixAlgebra().multiply(toIntegralBasisBaseChange(),
				field.asVector(t));
		List<IntE> rounded = new ArrayList<>();
		for (Fraction coefficient : inIntegralBasis.asList()) {
			rounded.add(coefficient.round());
		}
		return fromVector(new Vector<>(rounded));
	}

	public Matrix<Fraction> toIntegralBasisBaseChange() {
		getModuleGenerators();
		return toIntegralBasis;
	}

	@Override
	public Value valuation(NFE t, Ideal<NFE> maximalIdeal) {
		return maximalIdeal.maximumPowerContains(t);
	}

	@Override
	public NumberField fieldOfFractions() {
		return field;
	}

	@Override
	public MathMap<NFE, NFE> embedding() {
		return new Identity<>();
	}

	@Override
	public boolean isInteger(NFE t) {
		return isElement(t);
	}

	@Override
	public NFE asInteger(NFE t) {
		if (!isElement(t)) {
			throw new ArithmeticException("Not an integer!");
		}
		return t;
	}

	@Override
	public FiniteField reduction(Ideal<NFE> maximalIdeal) {
		localize(maximalIdeal);
		return localizations.get(maximalIdeal).reduction();
	}

	@Override
	public FFE reduce(NFE t, Ideal<NFE> maximalIdeal) {
		return localize(maximalIdeal).reduce(t);
	}

	@Override
	public NFE lift(FFE s, Ideal<NFE> maximalIdeal) {
		return localize(maximalIdeal).lift(s);
	}

	@Override
	public LocalRing<NFE, FFE> localize(Ideal<NFE> maximalIdeal) {
		if (!localizations.containsKey(maximalIdeal)) {
			NumberFieldIdeal ideal = (NumberFieldIdeal) maximalIdeal;
			localizations.put(maximalIdeal,
					new LocalizedNumberField(field, ideal, Integers.z().localize(ideal.prime().getValue())));
		}
		return localizations.get(maximalIdeal).ringOfIntegers();
	}

	@Override
	public NFE getEmbedding(IntE t) {
		return field.getEmbedding(t);
	}

	@Override
	public boolean isGeneratingAlgebra(List<NFE> s) {
		List<NFE> includingPowers = new ArrayList<>();
		includingPowers.add(one());
		for (NFE e : s) {
			NFE power = one();
			for (int i = 1; i < field.degree(); i++) {
				power = multiply(power, e);
				includingPowers.add(power);
			}
		}
		return isGeneratingModule(includingPowers);
	}

	@Override
	public List<NFE> getAlgebraGenerators() {
		return getModuleGenerators();
	}

	public class NumberFieldIdeal extends AbstractIdeal<NFE> implements Ideal<NFE> {
		private NFE uniformizer;
		private IntE intGenerator;
		private Matrix<IntE> inBasis;
		private FreeSubModule<IntE, NFE> asSubModule;
		private MatrixModule<IntE> matrixModule;
		private List<NFE> generators;
		private LocalRing<Fraction, PFE> localized;
		private OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type;
		private boolean maximal;
		private Map<IntE, FactorizationResult<Ideal<NFE>, Ideal<NFE>>> idealFactorsPerPrime;
		private Map<Ideal<NFE>, NFE> chineseRemainderMultipliers;
		private FactorizationResult<Ideal<NFE>, Ideal<NFE>> idealFactorization;
		private boolean zeroIdeal;
		private boolean unitIdeal;

		private NumberFieldIdeal(boolean zeroIdeal) {
			super(NumberFieldIntegers.this);
			if (zeroIdeal) {
				this.zeroIdeal = true;
				this.generators = Collections.emptyList();
				return;
			}
			this.unitIdeal = true;
			this.uniformizer = one();
			this.generators = Collections.singletonList(uniformizer);
			this.idealFactorization = new FactorizationResult<>(this, Collections.emptySortedMap());
		}

		private NumberFieldIdeal(LocalRing<Fraction, PFE> localized,
				List<OkutsuType<Fraction, PFE, PFE, FFE, FiniteField>> types, int index,
				List<List<UnivariatePolynomial<Fraction>>> integral) {
			super(NumberFieldIntegers.this);
			this.localized = localized;
			this.type = types.get(index);
			this.maximal = true;
			Rationals q = Rationals.q();
			UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
			UnivariatePolynomial<Fraction> uniformizer = null;
			UnivariatePolynomial<Fraction> basisElement = null;
			if (type.ramificationIndex() == 1) {
				uniformizer = polynomials.getEmbedding(localized.uniformizer());
				basisElement = integral.get(index).get(0);
			} else {
				boolean foundUniformizer = false;
				boolean foundZero = false;
				for (int i = 0; i < integral.get(index).size(); i++) {
					Value value = type.valuation(integral.get(index).get(i));
					if (!foundUniformizer && value.equals(Value.ONE)) {
						foundUniformizer = true;
						uniformizer = integral.get(index).get(i);
					} else if (!foundZero && value.equals(Value.ZERO)) {
						foundZero = true;
						basisElement = integral.get(index).get(i);
					}
					if (foundZero && foundUniformizer) {
						break;
					}
				}
				if (!foundUniformizer) {
					System.err.println("Did not find uniformizer, but kinda expected to find one!");
					uniformizer = type.lift(type.reduction().extension().one(), 1);
				}
				if (!foundZero) {
					throw new ArithmeticException("Did not find zero valuation element in integral basis!");
				}
			}
			UnivariatePolynomial<Fraction> generator = polynomials.multiply(uniformizer, basisElement);
			for (int i = 0; i < types.size(); i++) {
				if (i == index) {
					continue;
				}
				OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> otherType = types.get(i);
				for (int j = 0; j < integral.get(i).size(); j++) {
					if (otherType.valuation(integral.get(i).get(j)).equals(Value.ZERO)) {
						generator = polynomials.add(generator, integral.get(i).get(j));
						break;
					}
				}
			}
			this.uniformizer = field.fromPolynomial(generator);
			this.intGenerator = localized.uniformizer().asInteger();
			this.idealFactorization = new FactorizationResult<>(getUnitIdeal(), SingletonSortedMap.map(this, 1));
			init();
		}

		private NumberFieldIdeal(Map<IntE, FactorizationResult<Ideal<NFE>, Ideal<NFE>>> factorization) {
			super(NumberFieldIntegers.this);
			this.idealFactorsPerPrime = factorization;
			SortedMap<Ideal<NFE>, Integer> factors = new TreeMap<>();
			Map<IntE, Integer> hMap = new TreeMap<>();
			Map<IntE, NFE> alphaMap = new TreeMap<>();
			for (IntE prime : factorization.keySet()) {
				FactorizationResult<Ideal<NFE>, Ideal<NFE>> primeFactorization = factorization.get(prime);
				NFE alpha = one();
				int h = 0;
				for (Ideal<NFE> primeIdeal : primeFactorization.primeFactors()) {
					NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
					if (!prime.equals(ideal.intGenerator)) {
						throw new ArithmeticException("Prime does match prime ideal!");
					}
					if (!ideal.maximal) {
						throw new ArithmeticException("Prime Ideal factor not prime!");
					}
					int multiplicity = primeFactorization.multiplicity(primeIdeal);
					factors.put(primeIdeal, multiplicity);
					alpha = multiply(alpha, power(ideal.uniformizer, multiplicity));
					h = Math.max(h, MiscAlgorithms.DivRoundUp(multiplicity, ideal.type.ramificationIndex()));
				}
				hMap.put(prime, h);
				alphaMap.put(prime, alpha);
			}
			Integers z = Integers.z();
			this.intGenerator = z.one();
			this.uniformizer = zero();
			for (IntE prime : factorization.keySet()) {
				this.intGenerator = z.multiply(intGenerator, z.power(prime, hMap.get(prime)));
				NFE alpha = alphaMap.get(prime);
				for (IntE otherPrime : factorization.keySet()) {
					if (prime.equals(otherPrime)) {
						continue;
					}
					alpha = multiply(z.power(otherPrime, hMap.get(otherPrime) + 1), alpha);
				}
				this.uniformizer = add(this.uniformizer, alpha);
			}
			this.idealFactorization = new FactorizationResult<>(getUnitIdeal(), factors);
			this.maximal = false;
			if (idealFactorization.isIrreducible()) {
				throw new ArithmeticException("Prime Ideals need to use different constructor!");
			}
			init();

		}

		private void init() {
			Vector<IntE> asVector = NumberFieldIntegers.this.asVector(uniformizer);
			List<IntE> reduced = new ArrayList<>();
			for (IntE c : asVector.asList()) {
				reduced.add(new IntE(c.getValue().mod(intGenerator.getValue())));
			}
			this.uniformizer = NumberFieldIntegers.this.fromVector(new Vector<>(reduced));
			if (this.uniformizer.equals(zero())) {
				this.uniformizer = field.getEmbedding(intGenerator);
			}
			List<NFE> asList = new ArrayList<>();
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(b, uniformizer);
				asList.add(multiplied);
			}
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(intGenerator, b);
				asList.add(multiplied);
			}
			asSubModule = new FreeSubModule<>(NumberFieldIntegers.this, asList);
			if (!asSubModule.contains(uniformizer)) {
				throw new ArithmeticException("Uniformizer not in submodule?");
			}
			if (!asSubModule.contains(field.getEmbedding(intGenerator))) {
				throw new ArithmeticException("int generator not in submodule?");
			}
			boolean skipLLL = false;
			for (NFE b : asSubModule.getBasis()) {
				skipLLL = skipLLL || field.norm(b).getNumerator().getValue().bitLength() > 512;
			}
			List<NFE> basis = asSubModule.getBasis();
			if (!skipLLL) {
				basis = field.latticeReduction(basis);
			}
			asSubModule = new FreeSubModule<>(NumberFieldIntegers.this, basis);
			if (!asSubModule.contains(uniformizer)) {
				throw new ArithmeticException("Uniformizer not in submodule?");
			}
			if (!asSubModule.contains(field.getEmbedding(intGenerator))) {
				throw new ArithmeticException("int generator not in submodule?");
			}
			matrixModule = new MatrixModule<>(Integers.z(), field.degree(), 2 * field.degree());
			BigInteger norm = field.norm(uniformizer).asInteger().getValue().abs();
			for (NFE b : basis) {
				if (isElement(field.divide(b, field.getEmbedding(intGenerator)))) {
					continue;
				}
				BigInteger basisNorm = field.norm(b).asInteger().getValue().abs();
				if (norm == null || norm.compareTo(basisNorm) >= 0) {
//					System.err.println("LLL improved");
					uniformizer = b;
					norm = basisNorm;
				}
			}
			if (field.degree() == 2 && field.discriminant().compareTo(Integers.z().zero()) > 0) {
				IntE idealNorm = norm();
				if (!idealNorm.getValue().equals(norm)) {
					Integers z = Integers.z();
					Rationals q = Rationals.q();
					IntE discriminant = field.discriminant();
					boolean mod23;
					if (discriminant.getValue().mod(BigInteger.valueOf(4)).equals(BigInteger.ZERO)) {
						discriminant = Integers.z().divideChecked(discriminant, z.getInteger(4));
						mod23 = true;
					} else {
						mod23 = false;
					}
					List<Vector<IntE>> pells = MiscAlgorithms.pellsEquation(discriminant,
							mod23 ? idealNorm : z.multiply(4, idealNorm), true);
					pellLoop: for (Vector<IntE> pell : pells) {
						NFE alpha = field.sqrt(field.getInteger(discriminant)).keySet().iterator().next();
						int[][] pm1 = new int[][] { new int[] { 1, 1 }, new int[] { -1, 1 }, new int[] { 1, -1 } };
						Fraction x = q.getEmbedding(pell.get(1));
						Fraction y = q.getEmbedding(pell.get(2));
						if (!mod23) {
							x = q.divide(x, q.getInteger(2));
							y = q.divide(y, q.getInteger(2));
						}
						for (int[] i : pm1) {
							NFE candidate = field.add(field.multiply(i[0], field.getEmbedding(x)),
									field.multiply(i[1], field.getEmbedding(y), alpha));
							if (isElement(field.divide(this.uniformizer, candidate))
									&& isElement(field.divide(field.getInteger(this.intGenerator), candidate))) {
								this.uniformizer = candidate;
								break pellLoop;
							}
						}
					}
				}
			}
//			System.err.println("LLL done");
			if (isElement(field.divide(uniformizer, field.getEmbedding(intGenerator)))) {
				this.uniformizer = field.getEmbedding(intGenerator);
			}
			List<Vector<IntE>> asVectorList = new ArrayList<>();
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(b, uniformizer);
				asVectorList.add(NumberFieldIntegers.this.asVector(multiplied));
			}
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(intGenerator, b);
				asVectorList.add(NumberFieldIntegers.this.asVector(multiplied));
			}
			this.inBasis = Matrix.fromColumns(asVectorList);
			matrixModule = new MatrixModule<>(Integers.z(), field.degree(), 2 * field.degree());
			this.generators = new ArrayList<>();
			this.generators.add(uniformizer);
			if (!isElement(field.divide(field.getEmbedding(intGenerator), uniformizer))) {
				this.generators.add(field.getEmbedding(intGenerator));
			}
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
		public boolean isPrime() {
			return zeroIdeal || isMaximal();
		}

		@Override
		public boolean isMaximal() {
			return maximal;
		}

		public IntE intGenerator() {
			return intGenerator;
		}

		public NFE uniformizer() {
			return uniformizer;
		}

		@Override
		public List<NFE> generators() {
			return generators;
		}

		@Override
		public List<NFE> generate(NFE t) {
			if (zeroIdeal) {
				return Collections.emptyList();
			}
			NFE inIdeal = field.subtract(t, residue(t));
			if (generators.size() == 1) {
				return Collections.singletonList(field.divide(inIdeal, uniformizer));
			}
			Vector<IntE> inIntegerBasis = NumberFieldIntegers.this.asVector(inIdeal);
			Vector<IntE> inIdealBasis = matrixModule.solve(inBasis, inIntegerBasis);
			List<NFE> g = new ArrayList<>();
			g.add(NumberFieldIntegers.this.fromVector(new Vector<>(inIdealBasis.asList().subList(0, field.degree()))));
			g.add(NumberFieldIntegers.this
					.fromVector(new Vector<>(inIdealBasis.asList().subList(field.degree(), 2 * field.degree()))));
			return g;
		}

		public NFE round(NFE t, int accuracy) {
			Value neededValue = new Value(accuracy);
			NFE rounded = zero();
			while (true) {
				Value value = type.valuation(t.asPolynomial());
				if (value.compareTo(neededValue) >= 0) {
					return rounded;
				}
				FFE reduced = type.reduce(t.asPolynomial());
				NFE lifted = field.fromPolynomial(type.lift(reduced, value.value()));
				rounded = add(rounded, lifted);
				t = subtract(t, lifted);
			}
		}

		@Override
		public NFE residue(NFE t) {
			if (maximal) {
				Value value = type.valuation(t.asPolynomial());
				if (value.compareTo(Value.ZERO) > 0) {
					return zero();
				}
				FFE reduced = type.reduce(t.asPolynomial());
				return field.fromPolynomial(type.lift(reduced, 0));
			}
			if (zeroIdeal) {
				return t;
			}
			if (unitIdeal) {
				return zero();
			}
			Map<Ideal<NFE>, NFE> multipliers = chineseRemainderMultipliers();
			NFE result = zero();
			for (Ideal<NFE> ideal : idealFactorization.primeFactors()) {
				result = add(result, multiply(multipliers.get(ideal),
						((NumberFieldIdeal) ideal).round(t, idealFactorization.multiplicity(ideal))));
			}
			return result;
		}

		private Map<Ideal<NFE>, NFE> chineseRemainderMultipliers() {
			if (chineseRemainderMultipliers == null) {
				chineseRemainderMultipliers = new TreeMap<>();
				Integers z = Integers.z();
				Rationals q = Rationals.q();
				IntE allPrimes = z.one();
				Map<IntE, IntE> powerMap = new TreeMap<>();
				for (IntE prime : idealFactorsPerPrime.keySet()) {
					FactorizationResult<Ideal<NFE>, Ideal<NFE>> idealFactors = idealFactorsPerPrime.get(prime);
					int h = 0;
					for (Ideal<NFE> primeIdeal : idealFactors.primeFactors()) {
						NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
						int multiplicity = idealFactors.multiplicity(primeIdeal);
						h = Math.max(h, MiscAlgorithms.DivRoundUp(multiplicity, ideal.type.ramificationIndex()));
					}
					IntE primePower = z.power(prime, h);
					allPrimes = z.multiply(primePower, allPrimes);
					powerMap.put(prime, primePower);
				}
				for (IntE prime : idealFactorsPerPrime.keySet()) {
					IntE primePower = powerMap.get(prime);
					IntE allOther = z.divideChecked(allPrimes, primePower);
					IntE primeMultiplier = z.multiply(allOther.getValue().modInverse(primePower.getValue()), allOther);
					FactorizationResult<Ideal<NFE>, Ideal<NFE>> idealFactors = idealFactorsPerPrime.get(prime);
					for (Ideal<NFE> ideal : idealFactors.primeFactors()) {
						NumberFieldIdeal primeIdeal = (NumberFieldIdeal) ideal;
						int multiplicity = idealFactors.multiplicity(ideal);
						int ramification = primeIdeal.type.ramificationIndex();
						NFE multiplier = one();
						Fraction multiplierValue = q.zero();
						for (Ideal<NFE> otherIdeal : idealFactors.primeFactors()) {
							NumberFieldIdeal otherPrimeIdeal = (NumberFieldIdeal) otherIdeal;
							if (otherPrimeIdeal == primeIdeal) {
								continue;
							}
							int otherMultiplicity = idealFactors.multiplicity(otherIdeal);
							int otherRamificiation = otherPrimeIdeal.type.ramificationIndex();
							Fraction neededOtherValue = q.getFraction(otherMultiplicity, otherRamificiation);
							Value multiplierV = otherPrimeIdeal.type.valuation(multiplier.asPolynomial());
							if (multiplierV.isInfinite() || multiplierV.value() > otherMultiplicity) {
								continue;
							}
							Fraction otherMultiplierValue = q.getFraction(multiplierV.value(), otherRamificiation);
							Value otherV = otherPrimeIdeal.type.valuation(otherPrimeIdeal.type.representative());
							int value = primeIdeal.type.valuation(otherPrimeIdeal.type.representative()).value();
							multiplierValue = q.add(multiplierValue, q.getFraction(value, ramification));
							if (otherV.isInfinite()) {
								multiplier = multiply(field.fromPolynomial(otherPrimeIdeal.type.representative()),
										multiplier);
								continue;
							}
							Fraction otherValue = q.add(otherMultiplierValue,
									q.getFraction(otherV.value(), otherRamificiation));
							UnivariatePolynomial<Fraction> representative = otherPrimeIdeal.type.representative();
							if (otherValue.compareTo(
									q.add(neededOtherValue, q.getEmbedding(multiplierValue.roundUp()))) < 0) {
								int needed = q
										.multiply(otherRamificiation,
												q.add(neededOtherValue, q.getEmbedding(multiplierValue.roundUp())))
										.roundUp().intValueExact();
								representative = otherPrimeIdeal.localized
										.singleFactorLifting(otherPrimeIdeal.type, needed).representative();
							}
							multiplier = multiply(field.fromPolynomial(representative), multiplier);
						}
						int mv = primeIdeal.type.valuation(multiplier.asPolynomial()).value();
						int mod = mv % ramification;
						int power = mv / ramification;
						if (mod != 0) {
							for (NFE basis : getModuleGenerators()) {
								Value basisValue = primeIdeal.type.valuation(basis.asPolynomial());
								if (!basisValue.isInfinite() && basisValue.value() == ramification - mod) {
									multiplier = field.multiply(basis, multiplier);
									break;
								}
							}
							power += 1;
						}
						multiplier = field.multiply(field.getEmbedding(q.power(q.getEmbedding(prime), -power)),
								multiplier);
						multiplier = multiply(field.fromPolynomial(primeIdeal.localized
								.invertInType(multiplier.asPolynomial(), primeIdeal.type, multiplicity * ramification)),
								multiplier);
						multiplier = multiply(primeMultiplier, multiplier);
						if (!primeIdeal.residue(multiplier).equals(one())) {
							throw new ArithmeticException("Multiplier calculation went wrong!");
						}
						chineseRemainderMultipliers.put(primeIdeal, multiplier);
					}
				}
			}
			return chineseRemainderMultipliers;
		}

		@Override
		public boolean contains(NFE t) {
			if (t.equals(zero()) || unitIdeal) {
				return true;
			}
			if (zeroIdeal) {
				return false;
			}
			for (Ideal<NFE> primeIdeal : idealFactorization.primeFactors()) {
				NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
				Value value = ideal.type.valuation(t.asPolynomial());
				if (value.compareTo(new Value(idealFactorization.multiplicity(primeIdeal))) < 0) {
					return false;
				}
			}
			return true;
		}

		@Override
		public Value maximumPowerContains(NFE t) {
			if (t.equals(zero()) || unitIdeal) {
				return Value.INFINITY;
			}
			if (zeroIdeal) {
				return Value.ZERO;
			}
			Value value = Value.INFINITY;
			for (Ideal<NFE> primeIdeal : idealFactorization.primeFactors()) {
				NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
				value = value.min(new Value(
						ideal.type.valuation(t.asPolynomial()).value() / idealFactorization.multiplicity(primeIdeal)));
			}
			return value;
		}

		@Override
		public boolean contains(Ideal<NFE> other) {
			NumberFieldIdeal otherIdeal = (NumberFieldIdeal) other;
			if (otherIdeal.zeroIdeal) {
				return true;
			}
			if (zeroIdeal) {
				return false;
			}
			if (unitIdeal) {
				return true;
			}
			if (otherIdeal.unitIdeal) {
				return false;
			}
			if (otherIdeal.maximal && !maximal) {
				return false;
			}
			if (otherIdeal.maximal && otherIdeal.uniformizer.equals(uniformizer)
					&& otherIdeal.intGenerator.equals(intGenerator)) {
				return true;
			}
			if (otherIdeal.maximal) {
				return false;
			}
			for (Ideal<NFE> otherPrimeFactor : otherIdeal.idealFactorization.primeFactors()) {
				if (idealFactorization.multiplicity(otherPrimeFactor) < otherIdeal.idealFactorization
						.multiplicity(otherPrimeFactor)) {
					return false;
				}
			}
			return true;
		}

		@Override
		public int compareTo(Ideal<NFE> other) {
			NumberFieldIdeal otherIdeal = (NumberFieldIdeal) other;
			if (!maximal || !otherIdeal.maximal) {
				return super.compareTo(other);
			}
			int cmp = intGenerator.compareTo(otherIdeal.intGenerator);
			if (cmp != 0) {
				return cmp;
			}
			return uniformizer.compareTo(otherIdeal.uniformizer);
		}

		@Override
		public boolean equalsIdeal(Ideal<NFE> other) {
			NumberFieldIdeal otherIdeal = (NumberFieldIdeal) other;
			if (!maximal || !otherIdeal.maximal) {
				return super.equalsIdeal(other);
			}
			return compareTo(other) == 0;
		}

		public IntE prime() {
			if (!maximal) {
				throw new ArithmeticException("not a maximal ideal!");
			}
			return intGenerator;
		}

		public OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type() {
			if (!maximal) {
				throw new ArithmeticException("not a maximal ideal!");
			}
			return type;
		}

		public IntE norm() {
			Integers z = Integers.z();
			if (maximal) {
				return z.getInteger(type.reduction().extension().getNumberOfElements());
			}
			IntE norm = z.one();
			for (Ideal<NFE> primeIdeal : idealFactorization.primeFactors()) {
				NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
				norm = z.multiply(norm, z.power(ideal.norm(), idealFactorization.multiplicity(primeIdeal)));
			}
			return norm;
		}

		public IntE discriminant() {
			return NumberFieldIntegers.this.discriminant(asSubModule.getBasis());
		}
	}
}
