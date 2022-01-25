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
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractIdeal;
import fields.helper.FieldEmbedding;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Integers.IntegerIdeal;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Algebra;
import fields.interfaces.DedekindRing;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.DiscreteValuationRing.OkutsuType;
import fields.interfaces.DiscreteValuationRing.TheMontesResult;
import fields.interfaces.Ideal;
import fields.interfaces.Lattice;
import fields.interfaces.MathMap;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.Value;
import fields.numberfields.NumberField.NFE;
import fields.vectors.DualVectorSpace;
import fields.vectors.DualVectorSpace.Dual;
import fields.vectors.FreeModule;
import fields.vectors.FreeSubModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Polytope;
import fields.vectors.RealLattice;
import fields.vectors.Vector;
import util.Identity;
import util.MiscAlgorithms;
import util.SingletonSortedMap;

public class NumberFieldIntegers extends AbstractAlgebra<IntE, NFE>
		implements Algebra<IntE, NFE>, DedekindRing<NFE, NFE, FFE>, Lattice<NFE, Real, Vector<Real>> {
	private NumberField field;
	private UnivariatePolynomial<IntE> minimalPolynomial;
	private List<NFE> integralBasis;
	private Matrix<Fraction> toIntegralBasis;
	private Matrix<Real> generatorsAsMatrix;
	private Matrix<Fraction> fromIntegralBasis;
	private Set<NFE> units;
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
	public boolean isReduced() {
		return true;
	}

	@Override
	public boolean isIrreducible() {
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
		NFE residue = divisorIdeal.residue(dividend);
		return new QuotientAndRemainderResult<>(field.divide(subtract(dividend, residue), divisor), residue);
	}

	@Override
	public BigInteger euclidMeasure(NFE t) {
		throw new ArithmeticException("Not a Euclidean ring!");
	}

	@Override
	public NFE projectToUnit(NFE t) {
		NumberFieldIdeal ideal = getIdeal(Collections.singletonList(t));
		if (ideal.generators().size() != 1) {
			throw new ArithmeticException("Could not project to unit!");
		}
		return divideChecked(t, ideal.generators().get(0));
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
				Set<NFE> zeta3 = field.roots(eisenstein).keySet();
				units.addAll(zeta3);
				for (NFE unit : zeta3) {
					units.add(field.negative(unit));
				}
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

	public SortedMap<IntE, NFE> primitiveRootsOfUnity() {
		Integers z = Integers.z();
		UnivariatePolynomialRing<NFE> polynomials = field.getUnivariatePolynomialRing();
		SortedMap<IntE, NFE> result = new TreeMap<>();
		result.put(z.one(), one());
		List<List<IntE>> primePowers = new ArrayList<>();
		for (IntE prime : z.setOfPrimes()) {
			if (prime.compareTo(z.getInteger(field.degree() + 1)) > 0) {
				break;
			}
			UnivariatePolynomial<NFE> unity = polynomials.toUnivariate(polynomials.divideChecked(
					polynomials.subtract(polynomials.getVarPower(prime.intValueExact()), polynomials.one()),
					polynomials.subtract(polynomials.getVar(), polynomials.one())));
			Map<NFE, Integer> roots = field.roots(unity);
			if (roots.isEmpty()) {
				continue;
			}
			NFE rootOfUnity = roots.keySet().iterator().next();
			result.put(prime, rootOfUnity);
			List<IntE> powers = new ArrayList<>();
			powers.add(z.one());
			powers.add(prime);
			for (int i = 2; z.multiply(prime.intValueExact() - 1, z.power(prime, i - 1)).intValueExact() <= field
					.degree(); i++) {
				roots = field.roots(rootOfUnity, prime.intValueExact());
				if (roots.isEmpty()) {
					break;
				}
				IntE primePower = z.power(prime, i);
				rootOfUnity = roots.keySet().iterator().next();
				result.put(primePower, rootOfUnity);
				powers.add(primePower);
			}
			primePowers.add(powers);
		}
		List<List<IntE>> factors = MiscAlgorithms.crossProduct(primePowers);
		for (List<IntE> factorList : factors) {
			IntE factor = z.one();
			NFE rootOfUnity = one();
			for (IntE primePower : factorList) {
				factor = z.multiply(primePower, factor);
				rootOfUnity = multiply(result.get(primePower), rootOfUnity);
			}
			result.put(factor, rootOfUnity);
		}
		return result;
	}

	public Set<NFE> freeUnitGroupGenerators() {
//		int numberOfUnits = field.realEmbeddings().size() + field.complexEmbeddings().size() - 1;
//		if (numberOfUnits == 0) {
//			return Collections.emptySet();
//		}
//		if (field.degree() == 2) {
//			Integers z = Integers.z();
//			IntE discriminant = field.discriminant();
//			boolean mod23;
//			if (discriminant.getValue().mod(BigInteger.valueOf(4)).equals(BigInteger.ZERO)) {
//				discriminant = Integers.z().divideChecked(discriminant, z.getInteger(4));
//				mod23 = true;
//			} else {
//				mod23 = false;
//			}
//			List<Vector<IntE>> pells = MiscAlgorithms.pellsEquation(discriminant, z.getInteger(mod23 ? -1 : -4), false);
//			pells.addAll(MiscAlgorithms.pellsEquation(discriminant, z.getInteger(mod23 ? 1 : 4), false));
//			Vector<IntE> pell = pells.get(0);
//			NFE alpha = field.sqrt(field.getInteger(discriminant)).keySet().iterator().next();
//			NFE fundamental;
//			if (mod23) {
//				fundamental = add(getInteger(pell.get(1)), multiply(pell.get(2), alpha));
//			} else {
//				fundamental = divide(add(getInteger(pell.get(1)), multiply(pell.get(2), alpha)), getInteger(2));
//			}
//			return Collections.singleton(fundamental);
//		}
		if (units == null) {
			idealEquivalenceClass(getUnitIdeal());
		}
		return units;
	}

	public List<NFE> unitGroupGenerators() {
		List<NFE> result = new ArrayList<>();
		result.addAll(freeUnitGroupGenerators());
		SortedMap<IntE, NFE> rootsOfUnity = primitiveRootsOfUnity();
		result.add(rootsOfUnity.get(rootsOfUnity.lastKey()));
		return result;
	}

	@Override
	public int krullDimension() {
		return 1;
	}

	@Override
	public List<Ideal<NFE>> maximalPrimeIdealChain(Ideal<NFE> start) {
		if (start.equals(getZeroIdeal())) {
			return maximalPrimeIdealChain(start,
					idealsOver(Integers.z().getIdeal(Collections.singletonList(Integers.z().getInteger(2)))).get(0));
		}
		return Collections.singletonList(primaryDecomposition(start).getRadicals().get(0));
	}

	@Override
	public List<Ideal<NFE>> maximalPrimeIdealChain(Ideal<NFE> start, Ideal<NFE> end) {
		if (!end.contains(start) || !end.isPrime()) {
			throw new ArithmeticException("Invalid preconditions!");
		}
		if (start.equals(getZeroIdeal())) {
			if (end.equals(getZeroIdeal())) {
				return Collections.singletonList(getZeroIdeal());
			}
			List<Ideal<NFE>> result = new ArrayList<>();
			result.add(start);
			result.add(end);
			return result;
		}
		return Collections.singletonList(end);
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
		SortedMap<IntE, FactorizationResult<Ideal<NFE>, Ideal<NFE>>> factorizations = new TreeMap<>();
		SortedMap<Ideal<NFE>, Integer> allIdealFactors = new TreeMap<>();
		for (IntE prime : normFactorization.primeFactors()) {
			SortedMap<Ideal<NFE>, Integer> idealFactors = new TreeMap<>();
			List<NumberFieldIdeal> ideals = idealsOver(prime);
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
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = new FactorizationResult<>(getUnitIdeal(),
				allIdealFactors);
		if (!idealsByFactorization.containsKey(factorization)) {
			if (factorization.isIrreducible()) {
				idealsByFactorization.put(factorization, (NumberFieldIdeal) factorization.firstPrimeFactor());
			} else {
				idealsByFactorization.put(factorization, new NumberFieldIdeal(factorizations));
			}
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
	public NumberFieldIdeal intersect(Ideal<NFE> t1, Ideal<NFE> t2) {
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
		Map<NumberFieldIdeal, Integer> newFactors = new TreeMap<>();
		for (NumberFieldIdeal ideal : factors.keySet()) {
			int multiplicity = factors.get(ideal);
			if (multiplicity != 0) {
				newFactors.put(ideal, multiplicity);
			}
		}
		factors = newFactors;
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
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> asFactorizationResult = new FactorizationResult<>(getUnitIdeal(),
				cast);
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
	public DedekindRing<NFE, NFE, FFE> asDedekindRing() {
		return this;
	}

	@Override
	public PrimaryDecompositionResult<NFE, NumberFieldIdeal> primaryDecomposition(Ideal<NFE> t) {
		NumberFieldIdeal ideal = (NumberFieldIdeal) t;
		List<NumberFieldIdeal> primaries = new ArrayList<>();
		List<NumberFieldIdeal> radicals = new ArrayList<>();
		for (Ideal<NFE> prime : ideal.idealFactorization.primeFactors()) {
			NumberFieldIdeal primeIdeal = (NumberFieldIdeal) prime;
			primaries.add(power(prime, ideal.idealFactorization.multiplicity(prime)));
			radicals.add(primeIdeal);
		}
		return new PrimaryDecompositionResult<>(primaries, radicals);
	}

	@Override
	public ModuloMaximalIdealResult<NFE, FFE> moduloMaximalIdeal(Ideal<NFE> ideal) {
		return new ModuloMaximalIdealResult<>(this, ideal, reduction(ideal), new MathMap<>() {
			@Override
			public FFE evaluate(NFE t) {
				return reduce(t, ideal);
			}
		}, new MathMap<>() {
			@Override
			public NFE evaluate(FFE t) {
				return lift(t, ideal);
			}
		});
	}

	@Override
	public ModuloIdealResult<NFE, ?> moduloIdeal(Ideal<NFE> ideal) {
		if (ideal.isMaximal()) {
			ModuloMaximalIdealResult<NFE, FFE> result = moduloMaximalIdeal(ideal);
			return new ModuloIdealResult<>(this, ideal, result.getField(), result.getReduction(), result.getLift());
		}
		if (ideal.isPrime()) {
			return new ModuloIdealResult<>(this, ideal, this, new Identity<>(), new Identity<>());
		}
		ModuloNumberFieldIdeal modulo = new ModuloNumberFieldIdeal(this, (NumberFieldIdeal) ideal);
		return new ModuloIdealResult<>(this, ideal, modulo, new MathMap<>() {
			@Override
			public NFE evaluate(NFE t) {
				return modulo.getEmbedding(t);
			}

		}, new MathMap<>() {
			@Override
			public NFE evaluate(NFE t) {
				return t;
			}
		});
	}

	@Override
	public NumberFieldIdeal getUnitIdeal() {
		if (unitIdeal == null) {
			this.unitIdeal = new NumberFieldIdeal(false);
			this.idealsByFactorization.put(unitIdeal.idealFactorization, unitIdeal);
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

	@Override
	public NumberFieldIdeal getNilRadical() {
		return getZeroIdeal();
	}

	public List<NumberFieldIdeal> idealsOver(int prime) {
		return idealsOver(BigInteger.valueOf(prime));
	}

	public List<NumberFieldIdeal> idealsOver(BigInteger prime) {
		return idealsOver(Integers.z().getInteger(prime));
	}

	public List<NumberFieldIdeal> idealsOver(Ideal<IntE> integerIdeal) {
		return idealsOver(integerIdeal.generators().get(0));
	}

	public List<NumberFieldIdeal> idealsOver(IntE prime) {
		if (!Integers.z().isPrime(prime)) {
			throw new ArithmeticException("Prime number required!");
		}
		if (idealsOverPrime.containsKey(prime)) {
			return idealsOverPrime.get(prime);
		}
		DiscreteValuationRing<Fraction, PFE> zp = Integers.z().localize(prime);
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

	public List<NumberFieldIdeal> idealsOver(Ideal<NFE> lowerIdeal,
			FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		if (!fieldEmbedding.getField().equals(field)) {
			throw new ArithmeticException("incorrect field embedding");
		}
		NumberFieldIdeal lowerPrimeIdeal = (NumberFieldIdeal) lowerIdeal;
		NFE embeddedLowerUniformizer = fieldEmbedding.getEmbedding(lowerPrimeIdeal.uniformizer);
		Integers z = Integers.z();
		List<NumberFieldIdeal> overRationals = idealsOver(z.getIdeal(lowerPrimeIdeal.intGenerator));
		List<NumberFieldIdeal> result = new ArrayList<>();
		for (NumberFieldIdeal ideal : overRationals) {
			if (ideal.contains(embeddedLowerUniformizer)) {
				result.add(ideal);
			}
		}
		return result;
	}

	public IntegerIdeal intersectToIntegers(Ideal<NFE> ideal) {
		return Integers.z().getIdeal(Collections.singletonList(((NumberFieldIdeal) ideal).intGenerator));
	}

	public NumberFieldIdeal intersectToLowerField(Ideal<NFE> ideal,
			FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		NumberFieldIdeal numberFieldIdeal = (NumberFieldIdeal) ideal;
		NumberFieldIntegers lowerOrder = fieldEmbedding.getEmbeddedField().maximalOrder();
		List<NFE> lowerIntegralBasis = new ArrayList<>();
		for (NFE t : lowerOrder.getModuleGenerators()) {
			lowerIntegralBasis.add(fieldEmbedding.getEmbedding(t));
		}
		FreeSubModule<IntE, NFE> lowerIntegers = new FreeSubModule<>(this, lowerIntegralBasis);
		FreeSubModule<IntE, NFE> intersected = lowerIntegers.intersection(numberFieldIdeal.asSubModule);
		List<NFE> lowerGenerators = new ArrayList<>();
		for (NFE intersect : intersected.getBasis()) {
			Vector<NFE> asVector = fieldEmbedding.asVector(intersect);
			for (int i = 1; i < asVector.dimension(); i++) {
				if (!asVector.get(i + 1).equals(lowerOrder.zero())) {
					throw new ArithmeticException("Could not intersect");
				}
			}
			lowerGenerators.add(asVector.get(1));
		}
		return lowerOrder.getIdeal(lowerGenerators);
	}

	public NumberFieldIdeal extend(Ideal<IntE> ideal) {
		return getIdeal(Collections.singletonList(getEmbedding(ideal.generators().get(0))));
	}

	public NumberFieldIdeal extend(Ideal<NFE> ideal, FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		if (!fieldEmbedding.getField().equals(field)) {
			throw new ArithmeticException("Field Embedding incorrect");
		}
		return (NumberFieldIdeal) getIdealEmbedding(ideal, fieldEmbedding.getEmbeddingMap());
	}

	public IntegerIdeal norm(Ideal<NFE> ideal) {
		NumberFieldIdeal numberFieldIdeal = (NumberFieldIdeal) ideal;
		return Integers.z().getIdeal(Collections.singletonList(numberFieldIdeal.norm()));
	}

	public NumberFieldIdeal normOver(Ideal<NFE> ideal, FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		if (!fieldEmbedding.getField().equals(field)) {
			throw new ArithmeticException("Field Embedding incorrect");
		}
		NumberFieldIntegers lowerOrder = fieldEmbedding.getEmbeddedField().maximalOrder();
		if (ideal.equals(getZeroIdeal())) {
			return lowerOrder.getZeroIdeal();
		}
		if (ideal.equals(getUnitIdeal())) {
			return lowerOrder.getUnitIdeal();
		}
		NumberFieldIdeal numberFieldIdeal = (NumberFieldIdeal) ideal;
		if (numberFieldIdeal.isMaximal()) {
			NumberFieldIdeal intersected = intersectToLowerField(numberFieldIdeal, fieldEmbedding);
			int degree = numberFieldIdeal.type().residueDegree() / intersected.type().residueDegree();
			return lowerOrder.power(intersected, degree);
		}
		NumberFieldIdeal result = lowerOrder.getUnitIdeal();
		for (Ideal<NFE> primeIdeal : numberFieldIdeal.idealFactorization.primeFactors()) {
			result = lowerOrder.multiply(result, lowerOrder.power(normOver(primeIdeal, fieldEmbedding),
					numberFieldIdeal.idealFactorization.multiplicity(primeIdeal)));
		}
		return result;
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

	public NFE different(NFE t) {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> minimalPolynomial = field.minimalPolynomial(t);
		if (minimalPolynomial.degree() != field.degree()) {
			return zero();
		}
		UnivariatePolynomial<Fraction> derivative = polynomials.derivative(minimalPolynomial);
		UnivariatePolynomialRing<NFE> numberFieldPolynomials = field.getUnivariatePolynomialRing();
		UnivariatePolynomial<NFE> embeddedDerivative = numberFieldPolynomials.getEmbedding(derivative,
				field.getEmbeddingMap());
		return numberFieldPolynomials.evaluate(embeddedDerivative, t);
	}

	public NFE differentOver(NFE t, FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		UnivariatePolynomialRing<NFE> polynomials = fieldEmbedding.getEmbeddedField().getUnivariatePolynomialRing();
		UnivariatePolynomial<NFE> minimalPolynomial = fieldEmbedding.minimalPolynomial(t);
		if (minimalPolynomial.degree() != fieldEmbedding.relativeDegree()) {
			return zero();
		}
		UnivariatePolynomial<NFE> derivative = polynomials.derivative(minimalPolynomial);
		UnivariatePolynomialRing<NFE> numberFieldPolynomials = fieldEmbedding.getField().getUnivariatePolynomialRing();
		UnivariatePolynomial<NFE> embeddedDerivative = numberFieldPolynomials.getEmbedding(derivative,
				fieldEmbedding.getEmbeddingMap());
		return numberFieldPolynomials.evaluate(embeddedDerivative, t);
	}

	public NumberFieldIdeal different() {
		List<NFE> differentGenerator = new ArrayList<>();
		for (NFE generator : getModuleGenerators()) {
			for (int i = 1; i < field.degree(); i++) {
				differentGenerator.add(different(power(generator, i)));
			}
		}
		return getIdeal(differentGenerator);
	}

	public NumberFieldIdeal differentOver(FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		if (!field.equals(fieldEmbedding.getField())) {
			throw new ArithmeticException("needs to be called on upper order");
		}
		List<NFE> differentGenerator = new ArrayList<>();
		for (NFE generator : getModuleGenerators()) {
			for (int i = 1; i < field.degree(); i++) {
				differentGenerator.add(differentOver(power(generator, i), fieldEmbedding));
			}
		}
//				for (IntE prime : potentiallyRamifiedPrimes()) {
//			List<NumberFieldIdeal> ideals = idealsOver(Integers.z().getIdeal(Collections.singletonList(prime)));
//			for (NumberFieldIdeal ideal : ideals) {
//				differentGenerator.add(differentOver(ideal.uniformizer(), fieldEmbedding));
//			}
//		}
		return getIdeal(differentGenerator);
	}

	public IntE discriminant() {
		return discriminant(getModuleGenerators());
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

	public NumberFieldIdeal discriminantOver(FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		return differentOver(fieldEmbedding).normOver(fieldEmbedding);
	}

	private Set<IntE> potentiallyRamifiedPrimes() {
		Integers z = Integers.z();
		return z.uniqueFactorization(z.getUnivariatePolynomialRing().discriminant(minimalPolynomial())).primeFactors();
	}

	public IntE norm(NFE t) {
		return field.norm(t).asInteger();
	}

	public IntE trace(NFE t) {
		return field.trace(t).asInteger();
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
				DiscreteValuationRing<Fraction, PFE> zp = z.localize(prime.getValue());
				TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes = zp.theMontesAlgorithm(
						zpMinimalPolynomial, fp.getExtension(fp.getUnivariatePolynomialRing().getVar()));
				List<UnivariatePolynomial<Fraction>> integral = zp.triagonalizeIntegralBasis(zpMinimalPolynomial,
						zp.integralBasis(zpMinimalPolynomial, theMontes, false));
				List<Integer> valueList = new ArrayList<>();
				for (int i = 0; i < minimalPolynomial().degree(); i++) {
					valueList.add(-zp.localField().valuation(integral.get(i).leadingCoefficient()).value());
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
				if (!isElement(element)) {
					throw new ArithmeticException("Integral basis algorithm wrong, " + element + " is not integer");
				}
				integralBasis.add(element);
//				asVectors.add(field.asVector(element));
//				polynomialBasis.add(b);
			}
			this.integralBasis = sublatticeReduction(integralBasis);
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

	public Matrix<Fraction> toIntegralBasisBaseChange() {
		getModuleGenerators();
		return toIntegralBasis;
	}

	@Override
	public Value valuation(NFE t, Ideal<NFE> maximalIdeal) {
		return maximalIdeal.maximumPowerContains(t);
	}

	public NFE getNumerator(NFE t) {
		if (isElement(t)) {
			return t;
		}
		return field.multiply(t, getDenominator(t));
	}

	public NFE getDenominator(NFE t) {
		if (isElement(t)) {
			return one();
		}
		Integers z = Integers.z();
		Vector<Fraction> asVector = field.asVector(t);
		Matrix<Fraction> baseChange = toIntegralBasisBaseChange();
		Vector<Fraction> asIntegralVector = field.asVectorSpace().matrixAlgebra().multiply(baseChange, asVector);
		IntE denominator = z.one();
		for (Fraction c : asIntegralVector.asList()) {
			denominator = z.lcm(denominator, c.getDenominator());
		}
		return getEmbedding(denominator);
	}

	public NFE asInteger(NFE t) {
		if (isElement(t)) {
			return t;
		}
		throw new ArithmeticException("Not an integer!");
	}

	@Override
	public FieldOfFractionsResult<NFE, NFE> fieldOfFractions() {
		return new FieldOfFractionsResult<>(this, field, new Identity<>(), new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				return getNumerator(t);
			}
		}, new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				return getDenominator(t);
			}
		}, new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				return asInteger(t);
			}
		});
	}

	// TODO: Fix
	@Override
	public LocalizeResult<NFE, NFE, NFE, FFE> localizeAtIdeal(Ideal<NFE> primeIdeal) {
		localize(primeIdeal);
		LocalizedNumberField localized = localizations.get(primeIdeal);
		return new LocalizeResult<>(this, primeIdeal, localized.ringOfIntegers(), new Identity<>(), new MathMap<>() {
			@Override
			public NFE evaluate(NFE t) {
				Value v = localized.valuation(t);
				if (v.compareTo(Value.ZERO) >= 0) {
					return t;
				}
				return field.multiply(t, power(localized.uniformizer(), -v.value()));
			}
		}, new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				Value v = localized.valuation(t);
				if (v.compareTo(Value.ZERO) >= 0) {
					return one();
				}
				return power(localized.uniformizer(), -v.value());
			}
		}, new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				if (!isElement(t)) {
					throw new ArithmeticException("Not an integer!");
				}
				return t;
			}
		});
	}

	@Override
	public boolean isInteger(NFE t) {
		return isElement(t);
	}

	public NumberField numberField() {
		return field;
	}

	@Override
	public FiniteRealVectorSpace getVectorSpace() {
		return field.minkowskiEmbeddingSpace();
	}

	@Override
	public Vector<Real> embedding(NFE t) {
		return field.minkowskiEmbedding(t);
	}

	public List<NFE> sublatticeReduction(List<NFE> sublatticeBasis) {
		return getVectorSpace().latticeReduction(sublatticeBasis, this);
	}

	@Override
	public int rank() {
		return field.degree();
	}

	public NumberFieldOrder asOrder() {
		return field.getOrder(getModuleGenerators());
	}

	FractionalIdeal reducedRepresentative(Ideal<NFE> b) {
		IdealGroup idealGroup = field.idealGroup();
		FractionalIdeal ideal = idealGroup.getEmbedding(b);
		NFE minimal = field.minkowskiEmbeddingSpace().latticeReduction(ideal, 1.0).get(0);
		FractionalIdeal result = idealGroup.operate(ideal, idealGroup.getPrincipalIdeal(field.inverse(minimal)));
		return result;
	}

	List<NumberFieldIdeal> idealEquivalenceClass(Ideal<NFE> b) {
		boolean calculateUnits = units == null && field.realEmbeddings().size() + field.complexEmbeddings().size() >= 2;
		int k = 0;
		List<FractionalIdeal> idealList = new ArrayList<>();
		idealList.add(reducedRepresentative(b));
		List<NFE> multipliers = null;
		if (calculateUnits) {
			multipliers = new ArrayList<>();
			multipliers.add(one());
			units = new TreeSet<>();
		}
		IdealGroup idealGroup = field.idealGroup();
		while (k < idealList.size()) {
			for (NFE t : neighbors(idealList.get(k))) {
				FractionalIdeal nextIdeal = idealGroup.operate(idealList.get(k),
						idealGroup.getPrincipalIdeal(field.inverse(t)));
				boolean found = false;
				for (int j = 0; j < idealList.size(); j++) {
					if (nextIdeal.equals(idealList.get(j))) {
						if (calculateUnits) {
							NFE newUnit = field.multiply(t, field.divide(multipliers.get(k), multipliers.get(j)));
							if (!isRootOfUnity(newUnit)) {
								units.add(newUnit);
							}
						}
						found = true;
						break;
					}
				}
				if (!found) {
					idealList.add(nextIdeal);
					if (calculateUnits) {
						multipliers.add(field.multiply(t, multipliers.get(k)));
					}
				}
			}
			k++;
		}
		List<NumberFieldIdeal> ideals = new ArrayList<>();
		for (FractionalIdeal ideal : idealList) {
			ideals.add(ideal.clearDenominator().getFirst());
		}
		if (calculateUnits) {
			List<NFE> unitList = new ArrayList<>();
			unitList.addAll(units);
			units.clear();
			List<Vector<Real>> logarithmList = new ArrayList<>();
			for (NFE unit : unitList) {
				logarithmList.add(field.logRepresentation(unit));
			}
			System.out.println("logList: " + logarithmList);
			RealLattice lattice = new RealLattice(field.logRepresentationSpace(), logarithmList);
			List<Vector<IntE>> asIntVectors = new ArrayList<>();
			for (Vector<Real> logVector : logarithmList) {
				asIntVectors.add(lattice.asVector(logVector));
			}
			Matrix<IntE> asMatrix = Matrix.fromColumns(asIntVectors);
			MatrixModule<IntE> module = asMatrix.getModule(Integers.z());
			FreeModule<IntE> free = new FreeModule<>(Integers.z(), lattice.rank());
			for (int i = 0; i < lattice.rank(); i++) {
				Vector<IntE> exponents = module.solve(asMatrix, free.getUnitVector(i + 1));
				NFE unit = one();
				for (int j = 0; j < unitList.size(); j++) {
					unit = multiply(power(unitList.get(j), exponents.get(j + 1)), unit);
				}
				units.add(unit);
			}
		}
		return ideals;
	}

	private boolean isRootOfUnity(NFE unit) {
		Reals r = Reals.r(128);
		for (Real value : values(unit).asList()) {
			if (!r.close(value, r.one())) {
				return false;
			}
		}
		return true;
	}

	private Set<NFE> neighbors(FractionalIdeal ideal) {
		int l = 0;
		Reals r = Reals.r(128);
		List<Set<NFE>> minimalSets = new ArrayList<>();
		Set<NFE> startingSet = new TreeSet<>();
		List<Real> border = new ArrayList<>();
		List<Dual<Real, Vector<Real>>> constraints = new ArrayList<>();
		DualVectorSpace<Real, Vector<Real>> dualSpace = field.minkowskiEmbeddingSpace().getDual();
		for (int i = 0; i < field.degree(); i++) {
			border.add(r.one());
			border.add(r.one());
			constraints.add(dualSpace.getUnitVector(i + 1));
			constraints.add(dualSpace.negative(dualSpace.getUnitVector(i + 1)));
		}
		Real epsilon = r.getPowerOfTwo(-r.precision()/2);
		Real boundary = r.add(r.one(), epsilon);
		candidateLoop: for (NFE candidate : field.minkowskiEmbeddingSpace()
				.latticePointsInPolytope(new Polytope<>(field.minkowskiEmbeddingSpace(), constraints, border), ideal)) {
			if (candidate.equals(zero())) {
				continue;
			}
			Vector<Real> values = values(candidate);
			for (int i = 0; i < values.dimension(); i++) {
				if (values.get(i + 1).compareTo(boundary) > 0) {
					continue candidateLoop;
				}
			}
			startingSet.add(candidate);
		}
//		System.err.println("Ideal: " + ideal);
//		System.err.println("Starting Set: " + startingSet);
		minimalSets.add(startingSet);
		int embeddings = field.logRepresentationSpace().dimension();
		while (l < minimalSets.size()) {
//			System.out.println("l: " + l + "\nMinimal Sets: " + minimalSets);
			for (int i = 0; i < embeddings; i++) {
				Set<NFE> set = minimalSets.get(l);
				Set<NFE> expansion = expansion(set, i + 1, ideal);
//				System.out.println("Expansion " + (i + 1) + ": " + expansion);
				if (!minimalSets.contains(expansion)) {
//					System.err.println("Added Expansion: " + expansion);
					minimalSets.add(expansion);
				}
				if (value(set, i + 1).compareTo(r.one()) > 0) {
					Set<NFE> compression = compression(set, i + 1);
//					System.out.println("Compression " + (i + 1) + ": " + compression);
					if (!minimalSets.contains(compression)) {
//						System.err.println("Added Compression: " + compression);
						minimalSets.add(compression);
					}
				}
			}
			l++;
		}
		Set<NFE> result = new TreeSet<>();
		for (Set<NFE> minimalSet : minimalSets) {
			result.addAll(minimalSet);
		}
		return result;
	}

//	private Real logValue(Set<NFE> minimalSet, int index) {
//		Real result = null;
//		for (NFE alpha : minimalSet) {
//			Real value = field.logRepresentation(alpha).get(index);
//			if (result == null || value.compareTo(result) > 0) {
//				result = value;
//			}
//		}
//		return result;
//	}
//
//	private Vector<Real> logValues(Set<NFE> minimalSet) {
//		List<Real> result = new ArrayList<>();
//		for (int i = 0; i < field.realEmbeddings().size() + field.complexEmbeddings().size(); i++) {
//			result.add(logValue(minimalSet, i + 1));
//		}
//		return new Vector<>(result);
//	}

	private Real realValue(Set<NFE> minimalSet, EmbeddedNumberField<Real> embeddedField) {
		Real result = null;
		for (NFE alpha : minimalSet) {
			Real value = embeddedField.value(alpha);
			if (result == null || value.compareTo(result) > 0) {
				result = value;
			}
		}
		return result;
	}

	private Real complexValue(Set<NFE> minimalSet, EmbeddedNumberField<ComplexNumber> embeddedField) {
		Real result = null;
		Complex c = Complex.c(128);
		for (NFE alpha : minimalSet) {
			Real value = c.norm(embeddedField.embedding(alpha));
			if (result == null || value.compareTo(result) > 0) {
				result = value;
			}
		}
		return result;
	}

	private Real value(Set<NFE> minimalSet, int index) {
		if (index <= field.realEmbeddings().size()) {
			return realValue(minimalSet, field.realEmbeddings().get(index - 1));
		}
		return complexValue(minimalSet, field.complexEmbeddings().get(index - field.realEmbeddings().size() - 1));
	}

	private Vector<Real> values(Set<NFE> minimalSet) {
		List<Real> result = new ArrayList<>();
		for (int i = 0; i < field.realEmbeddings().size() + field.complexEmbeddings().size(); i++) {
			result.add(value(minimalSet, i + 1));
		}
		return new Vector<>(result);
	}

	private Real value(NFE alpha, int index) {
		if (index <= field.realEmbeddings().size()) {
			return field.realEmbeddings().get(index - 1).value(alpha);
		}
		EmbeddedNumberField<ComplexNumber> embedding = field.complexEmbeddings()
				.get(index - field.realEmbeddings().size() - 1);
		Complex c = Complex.c(128);
		return c.norm(embedding.embedding(alpha));
	}

	private Vector<Real> values(NFE alpha) {
		List<Real> result = new ArrayList<>();
		for (int i = 0; i < field.realEmbeddings().size() + field.complexEmbeddings().size(); i++) {
			result.add(value(alpha, i + 1));
		}
		return new Vector<>(result);
	}

	private Set<NFE> restriction(Set<NFE> minimalSet, int index) {
		Set<NFE> result = new TreeSet<>();
		Vector<Real> values = values(minimalSet);
		alphaLoop: for (NFE alpha : minimalSet) {
			Vector<Real> v = values(alpha);
			if (!v.get(index).equals(values.get(index))) {
				continue;
			}
			for (int i = 0; i < field.realEmbeddings().size() + field.complexEmbeddings().size(); i++) {
				if (i + 1 == index) {
					continue;
				}
				if (v.get(i + 1).compareTo(values.get(i + 1)) >= 0) {
					continue alphaLoop;
				}
			}
			result.add(alpha);
		}
		return result;
	}

	private Set<NFE> expansion(Set<NFE> minimalSet, int index, FractionalIdeal ideal) {
		if (!restriction(minimalSet, index).isEmpty()) {
			return minimalSet;
		}
		Reals r = Reals.r(128);
		Real epsilon = r.power(r.getInteger(2), -16);
		Vector<Real> absoluteValues = values(minimalSet);
		int realEmbeddings = field.realEmbeddings().size();
		int complexEmbeddings = field.complexEmbeddings().size();
		List<Real> values = new ArrayList<>();
		for (int i = 0; i < realEmbeddings + complexEmbeddings; i++) {
			if (i + 1 == index) {
				Real limit = r.abs(r.multiply(r.getInteger(100 * field.degree()), r.getInteger(field.discriminant())));
				values.add(limit);
				values.add(r.negative(epsilon));
				if (i >= realEmbeddings) {
					values.add(limit);
					values.add(r.negative(epsilon));
				}
			} else {
				Real value = absoluteValues.get(i + 1);
				if (i >= realEmbeddings) {
					value = r.multiply(r.inverseSqrt(r.getInteger(2)), r.positiveSqrt(value));
				}
				value = r.subtract(value, epsilon);
				values.add(value);
				values.add(value);
				if (i >= realEmbeddings) {
					values.add(value);
					values.add(value);
				}
			}
		}
		List<Dual<Real, Vector<Real>>> constraints = new ArrayList<>();
		DualVectorSpace<Real, Vector<Real>> dualVectorSpace = field.minkowskiEmbeddingSpace().getDual();
		for (int i = 0; i < field.degree(); i++) {
			constraints.add(dualVectorSpace.getUnitVector(i + 1));
			constraints.add(dualVectorSpace.negative(dualVectorSpace.getUnitVector(i + 1)));
		}
		Polytope<Real, Vector<Real>> polytope = new Polytope<>(field.minkowskiEmbeddingSpace(), constraints, values);
		int realIndex = index <= realEmbeddings ? index : realEmbeddings + 2 * (index - realEmbeddings) - 1;
		Dual<Real, Vector<Real>> minimize = index <= realEmbeddings
				? dualVectorSpace.negative(dualVectorSpace.getUnitVector(realIndex))
				: dualVectorSpace.add(dualVectorSpace.negative(dualVectorSpace.getUnitVector(realIndex)),
						dualVectorSpace.negative(dualVectorSpace.getUnitVector(realIndex + 1)));
		NFE betterBorder = field.minkowskiEmbeddingSpace().latticeLinearProgram(polytope, minimize, ideal);
		values.clear();
		for (int i = 0; i < realEmbeddings + complexEmbeddings; i++) {
			if (i + 1 == index) {
				Real limit = values(betterBorder).get(index);
				if (i >= realEmbeddings) {
					limit = r.positiveSqrt(limit);
				}
				values.add(limit);
				values.add(limit);
				if (i >= realEmbeddings) {
					values.add(limit);
					values.add(limit);
				}
			} else {
				Real value = absoluteValues.get(i + 1);
				if (i >= realEmbeddings) {
					value = r.positiveSqrt(value);
				}
				values.add(value);
				values.add(value);
				if (i >= realEmbeddings) {
					values.add(value);
					values.add(value);
				}
			}
		}
		polytope = new Polytope<>(field.minkowskiEmbeddingSpace(), constraints, values);
		List<NFE> candidates = field.minkowskiEmbeddingSpace().latticePointsInPolytope(polytope, ideal);
		Real minValue = null;
		Set<NFE> add = new TreeSet<>();
		epsilon = r.getPowerOfTwo(-r.precision()/2);
		candidateLoop: for (NFE candidate : candidates) {
			if (candidate.equals(zero())) {
				continue;
			}
			Vector<Real> candidateValues = values(candidate);
			for (int i = 0; i < realEmbeddings + complexEmbeddings; i++) {
				if (i + 1 == index) {
					continue;
				}
				if (candidateValues.get(i + 1).compareTo(r.subtract(absoluteValues.get(i + 1), epsilon)) >= 0) {
					continue candidateLoop;
				}
			}
			if (minValue == null || minValue.compareTo(candidateValues.get(index)) > 0) {
				add.clear();
				add.add(candidate);
				minValue = candidateValues.get(index);
			} else if (minValue.equals(candidateValues.get(index))) {
				add.add(candidate);
			}
		}
		Set<NFE> result = new TreeSet<>();
		result.addAll(minimalSet);
		result.addAll(add);
		return result;
	}

	private Set<NFE> compression(Set<NFE> minimalSet, int index) {
		Set<NFE> result = new TreeSet<>();
		Real value = value(minimalSet, index);
		for (NFE alpha : minimalSet) {
			Real v = value(alpha, index);
			if (v.compareTo(value) < 0) {
				result.add(alpha);
			}
		}
		return result;
	}

	@Override
	public FiniteField reduction(Ideal<NFE> maximalIdeal) {
		localize(maximalIdeal);
		return localizations.get(maximalIdeal).residueField();
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
	public DiscreteValuationRing<NFE, FFE> localize(Ideal<NFE> maximalIdeal) {
		if (!localizations.containsKey(maximalIdeal)) {
			NumberFieldIdeal ideal = (NumberFieldIdeal) maximalIdeal;
			localizations.put(maximalIdeal,
					new LocalizedNumberField(field, ideal, Integers.z().localize(ideal.prime().getValue())));
		}
		return localizations.get(maximalIdeal).ringOfIntegers();
	}

	public LocalizedNumberField localizeAndQuotient(Ideal<NFE> maximalIdeal) {
		localize(maximalIdeal);
		return localizations.get(maximalIdeal);
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
		private DiscreteValuationRing<Fraction, PFE> localized;
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
				this.uniformizer = zero();
				this.intGenerator = Integers.z().zero();
				this.generators = Collections.emptyList();
				return;
			}
			this.unitIdeal = true;
			this.uniformizer = one();
			this.intGenerator = Integers.z().one();
			this.generators = Collections.singletonList(uniformizer);
			this.idealFactorization = new FactorizationResult<>(this, Collections.emptySortedMap());
		}

		private NumberFieldIdeal(DiscreteValuationRing<Fraction, PFE> localized,
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
			for (int i = 0; i < types.size(); i++) {
				if (i == index) {
					if (!types.get(i).valuation(generator).equals(Value.ONE)) {
						throw new ArithmeticException("Not a uniformizer of the prime ideal!");
					}
				} else {
					if (!types.get(i).valuation(generator).equals(Value.ZERO)) {
						throw new ArithmeticException("Not a unit mod the other prime ideals!");
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

		private boolean proposeNewUniformizer(NFE proposedUniformizer) {
			if (!asSubModule.contains(proposedUniformizer)) {
				return false;
			}
			if (proposedUniformizer.asPolynomial().degree() <= 0) {
				return false;
			}
			if (field.norm(proposedUniformizer).asInteger().getValue().abs()
					.compareTo(field.norm(uniformizer).asInteger().getValue().abs()) >= 0) {
				return false;
			}
			List<NFE> asList = new ArrayList<>();
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(b, proposedUniformizer);
				asList.add(multiplied);
			}
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(intGenerator, b);
				asList.add(multiplied);
			}
			FreeSubModule<IntE, NFE> proposedAsSubModule = new FreeSubModule<>(NumberFieldIntegers.this, asList);
			if (!proposedAsSubModule.contains(uniformizer)) {
				return false;
			}
			this.uniformizer = proposedUniformizer;
			this.asSubModule = proposedAsSubModule;
			return true;
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
			matrixModule = new MatrixModule<>(Integers.z(), field.degree(), 2 * field.degree());
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
				basis = sublatticeReduction(basis);
			}
//			asSubModule = new FreeSubModule<>(NumberFieldIntegers.this, basis);
//			if (!asSubModule.contains(uniformizer)) {
//				throw new ArithmeticException("Uniformizer not in submodule?");
//			}
//			if (!asSubModule.contains(field.getEmbedding(intGenerator))) {
//				throw new ArithmeticException("int generator not in submodule?");
//			}
			BigInteger norm = field.norm(uniformizer).asInteger().getValue().abs();
			for (NFE b : basis) {
				proposeNewUniformizer(b);
//				if (isElement(field.divide(b, field.getEmbedding(intGenerator)))) {
//					continue;
//				}
//				BigInteger basisNorm = field.norm(b).asInteger().getValue().abs();
//				if (norm == null || norm.compareTo(basisNorm) >= 0) {
////					System.err.println("LLL improved");
//					uniformizer = b;
//					norm = basisNorm;
//				}
			}
//			if (field.degree() == 2 && field.discriminant().compareTo(Integers.z().zero()) > 0) {
//				IntE idealNorm = norm();
//				if (!idealNorm.getValue().equals(norm)) {
//					Integers z = Integers.z();
//					Rationals q = Rationals.q();
//					IntE discriminant = field.discriminant();
//					boolean mod23;
//					if (discriminant.getValue().mod(BigInteger.valueOf(4)).equals(BigInteger.ZERO)) {
//						discriminant = Integers.z().divideChecked(discriminant, z.getInteger(4));
//						mod23 = true;
//					} else {
//						mod23 = false;
//					}
//					List<Vector<IntE>> pells = MiscAlgorithms.pellsEquation(discriminant,
//							mod23 ? idealNorm : z.multiply(4, idealNorm), true);
//					pellLoop: for (Vector<IntE> pell : pells) {
//						NFE alpha = field.sqrt(field.getInteger(discriminant)).keySet().iterator().next();
//						int[][] pm1 = new int[][] { new int[] { 1, 1 }, new int[] { -1, 1 }, new int[] { 1, -1 } };
//						Fraction x = q.getEmbedding(pell.get(1));
//						Fraction y = q.getEmbedding(pell.get(2));
//						if (!mod23) {
//							x = q.divide(x, q.getInteger(2));
//							y = q.divide(y, q.getInteger(2));
//						}
//						for (int[] i : pm1) {
//							NFE candidate = field.add(field.multiply(i[0], field.getEmbedding(x)),
//									field.multiply(i[1], field.getEmbedding(y), alpha));
//							if (isElement(field.divide(this.uniformizer, candidate))
//									&& isElement(field.divide(field.getInteger(this.intGenerator), candidate))) {
//								proposeNewUniformizer(candidate);
//								break pellLoop;
//							}
//						}
//					}
//				}
//			}
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
		public boolean isPrincipal() {
			if (generators.size() == 1) {
				return true;
			}
			NumberFieldIdeal representative = reducedRepresentative(this).clearDenominator().getFirst();
			for (NumberFieldIdeal ideal : field.idealClassGroup().neutral().alternativeRepresentatives()) {
				if (ideal.equals(representative)) {
					return true;
				}
			}
			return false;
		}

		@Override
		public boolean isPrimary() {
			return isPrime() || idealFactorization.primeFactors().size() == 1;
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
		public List<List<NFE>> nonTrivialCombinations(List<NFE> s) {
			List<Vector<IntE>> asVectorList = new ArrayList<>();
			for (NFE e : s) {
				for (NFE basisVector : NumberFieldIntegers.this.getModuleGenerators()) {
					asVectorList.add(NumberFieldIntegers.this.asVector(multiply(e, basisVector)));
				}
			}
			List<List<IntE>> nonTrivialCombinations = asFreeModule.nonTrivialCombinations(asVectorList);
			List<List<NFE>> result = new ArrayList<>();
			int degree = field.degree();
			for (List<IntE> combination : nonTrivialCombinations) {
				List<NFE> row = new ArrayList<>();
				for (int i = 0; i < s.size(); i++) {
					row.add(NumberFieldIntegers.this
							.fromVector(new Vector<>(combination.subList(i * degree, (i + 1) * degree))));
				}
				boolean found = false;
				for (List<NFE> prevRow : result) {
					NFE ratio = null;
					boolean notFound = false;
					for (int i = 0; i < s.size(); i++) {
						QuotientAndRemainderResult<NFE> qr = quotientAndRemainder(row.get(i), prevRow.get(i));
						if (!qr.getRemainder().equals(zero())) {
							notFound = true;
							break;
						}
						if (ratio == null) {
							ratio = qr.getQuotient();
						}
						if (!ratio.equals(qr.getQuotient())) {
							notFound = true;
							break;
						}
					}
					found = !notFound;
					if (found) {
						break;
					}
				}
				if (!found) {
					result.add(row);
				}
			}
			return result;
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
				PrimaryDecompositionResult<NFE, NumberFieldIdeal> primaryDecomposition = primaryDecomposition(this);
				ChineseRemainderPreparation<NFE> crp = prepareChineseRemainderTheorem(
						primaryDecomposition.getPrimaryIdeals());
				for (int i = 0; i < primaryDecomposition.getRadicals().size(); i++) {
					Ideal<NFE> ideal = primaryDecomposition.getRadicals().get(i);
					chineseRemainderMultipliers.put(ideal, crp.getMultipliers().get(i));
				}
				return chineseRemainderMultipliers;
//				Integers z = Integers.z();
//				Rationals q = Rationals.q();
//				IntE allPrimes = z.one();
//				Map<IntE, IntE> powerMap = new TreeMap<>();
//				for (IntE prime : idealFactorsPerPrime.keySet()) {
//					FactorizationResult<Ideal<NFE>, Ideal<NFE>> idealFactors = idealFactorsPerPrime.get(prime);
//					int h = 0;
//					for (Ideal<NFE> primeIdeal : idealFactors.primeFactors()) {
//						NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
//						int multiplicity = idealFactors.multiplicity(primeIdeal);
//						h = Math.max(h, MiscAlgorithms.DivRoundUp(multiplicity, ideal.type.ramificationIndex()));
//					}
//					IntE primePower = z.power(prime, h);
//					allPrimes = z.multiply(primePower, allPrimes);
//					powerMap.put(prime, primePower);
//				}
//				for (IntE prime : idealFactorsPerPrime.keySet()) {
//					IntE primePower = powerMap.get(prime);
//					IntE allOther = z.divideChecked(allPrimes, primePower);
//					IntE primeMultiplier = z.multiply(allOther.getValue().modInverse(primePower.getValue()), allOther);
//					FactorizationResult<Ideal<NFE>, Ideal<NFE>> idealFactors = idealFactorsPerPrime.get(prime);
//					for (Ideal<NFE> ideal : idealFactors.primeFactors()) {
//						NumberFieldIdeal primeIdeal = (NumberFieldIdeal) ideal;
//						int multiplicity = idealFactors.multiplicity(ideal);
//						int ramification = primeIdeal.type.ramificationIndex();
//						NFE multiplier = one();
//						Fraction multiplierValue = q.zero();
//						for (Ideal<NFE> otherIdeal : idealsOver(z.getIdeal(Collections.singletonList(prime)))/*idealFactors.primeFactors()*/) {
//							NumberFieldIdeal otherPrimeIdeal = (NumberFieldIdeal) otherIdeal;
//							if (otherPrimeIdeal == primeIdeal) {
//								continue;
//							}
//							int otherMultiplicity = idealFactors.multiplicity(otherIdeal);
//							int otherRamificiation = otherPrimeIdeal.type.ramificationIndex();
//							Fraction neededOtherValue = q.getFraction(otherMultiplicity, otherRamificiation);
//							Value multiplierV = otherPrimeIdeal.type.valuation(multiplier.asPolynomial());
//							if (multiplierV.isInfinite() || multiplierV.value() > otherMultiplicity) {
//								continue;
//							}
//							Fraction otherMultiplierValue = q.getFraction(multiplierV.value(), otherRamificiation);
//							Value otherV = otherPrimeIdeal.type.valuation(otherPrimeIdeal.type.representative());
//							int value = primeIdeal.type.valuation(otherPrimeIdeal.type.representative()).value();
//							multiplierValue = q.add(multiplierValue, q.getFraction(value, ramification));
//							if (otherV.isInfinite()) {
//								multiplier = multiply(field.fromPolynomial(otherPrimeIdeal.type.representative()),
//										multiplier);
//								continue;
//							}
//							Fraction otherValue = q.add(otherMultiplierValue,
//									q.getFraction(otherV.value(), otherRamificiation));
//							UnivariatePolynomial<Fraction> representative = otherPrimeIdeal.type.representative();
//							if (otherValue.compareTo(
//									q.add(neededOtherValue, q.getEmbedding(multiplierValue.roundUp()))) < 0) {
//								int needed = q
//										.multiply(otherRamificiation,
//												q.add(neededOtherValue, q.getEmbedding(multiplierValue.roundUp())))
//										.roundUp().intValueExact();
//								representative = otherPrimeIdeal.localized
//										.singleFactorLifting(otherPrimeIdeal.type, needed).representative();
//							}
//							multiplier = multiply(field.fromPolynomial(representative), multiplier);
//						}
//						int mv = primeIdeal.type.valuation(multiplier.asPolynomial()).value();
//						int mod = mv % ramification;
//						int power = mv / ramification;
//						if (mod != 0) {
//							for (NFE basis : getModuleGenerators()) {
//								Value basisValue = primeIdeal.type.valuation(basis.asPolynomial());
//								if (!basisValue.isInfinite() && basisValue.value() == ramification - mod) {
//									multiplier = field.multiply(basis, multiplier);
//									break;
//								}
//							}
//							power += 1;
//						}
//						multiplier = field.multiply(field.getEmbedding(q.power(q.getEmbedding(prime), -power)),
//								multiplier);
//						multiplier = multiply(field.fromPolynomial(primeIdeal.localized
//								.invertInType(multiplier.asPolynomial(), primeIdeal.type, multiplicity * ramification)),
//								multiplier);
//						multiplier = multiply(primeMultiplier, multiplier);
//						if (!primeIdeal.residue(multiplier).equals(one())) {
//							throw new ArithmeticException(
//									"Multiplier calculation went wrong (does not reduce to one)!");
//						}
//						if (!primeIdeal.type.valuation(multiplier.asPolynomial()).equals(Value.ZERO)) {
//							throw new ArithmeticException("Multiplier calculation went wrong (Not zero value)!");
//						}
//						for (Ideal<NFE> otherIdeal : idealFactors.primeFactors()) {
//							NumberFieldIdeal otherPrimeIdeal = (NumberFieldIdeal) otherIdeal;
//							if (otherPrimeIdeal == primeIdeal) {
//								continue;
//							}
//							if (primeIdeal.type.valuation(multiplier.asPolynomial()).compareTo(Value.ZERO) <= 0) {
//								throw new ArithmeticException(
//										"Multiplier calculation went wrong (Not positive value)!");
//							}
//						}
//						chineseRemainderMultipliers.put(primeIdeal, multiplier);
//					}
//				}
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
				Set<Ideal<NFE>> ideals = new TreeSet<>();
				FactorizationResult<Ideal<NFE>, Ideal<NFE>> thisFactorization = idealFactorization(this);
				FactorizationResult<Ideal<NFE>, Ideal<NFE>> otherFactorization = idealFactorization(other);
				ideals.addAll(thisFactorization.primeFactors());
				ideals.addAll(otherFactorization.primeFactors());
				for (Ideal<NFE> ideal : ideals) {
					int cmp = thisFactorization.multiplicity(ideal) - otherFactorization.multiplicity(ideal);
					if (cmp != 0) {
						return cmp;
					}
				}
				return 0;
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

		public NumberFieldIdeal normOver(FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
			SortedMap<NumberFieldIdeal, Integer> factorization = new TreeMap<>();
			for (Ideal<NFE> primeFactor : idealFactorization(this).primeFactors()) {
				NumberFieldIdeal lowerPrimeIdeal = intersectToLowerField(primeFactor, fieldEmbedding);
				factorization.put(lowerPrimeIdeal, ((NumberFieldIdeal) primeFactor).type().residueDegree()
						/ lowerPrimeIdeal.type().residueDegree() * idealFactorization(this).multiplicity(primeFactor));
			}
			return fieldEmbedding.getEmbeddedField().maximalOrder().fromFactorization(factorization);
		}

		public IntE discriminant() {
			return NumberFieldIntegers.this.discriminant(asSubModule.getBasis());
		}

		public Real minkowskiBound() {
			Reals r = Reals.r(1024);
			Real twoByPiToS = r.power(r.divide(r.getInteger(2), r.pi()), field.complexEmbeddings().size());
			Real sqrtDiscriminant = r.positiveSqrt(r.abs(r.getEmbedding(field.discriminant())));
			Real norm = r.getEmbedding(norm());
			return r.multiply(twoByPiToS, sqrtDiscriminant, norm);
		}
	}
}
