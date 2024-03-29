package fields.finitefields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField.PFE;
import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.helper.FieldAutomorphism;
import fields.helper.FieldEmbedding;
import fields.helper.GaloisGroup;
import fields.helper.GenericAlgebraicRingExtension;
import fields.helper.GenericAlgebraicRingExtension.GenericAlgebraicExtensionElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.FieldExtension;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.Pair;

public class FiniteField extends AbstractFieldExtension<PFE, FFE, FiniteField>
		implements FieldExtension<PFE, FFE, FiniteField> {
	private static Map<BigInteger, FiniteField> finiteFields = new TreeMap<>();

	private Map<Polynomial<FFE>, List<Pair<Polynomial<FFE>, Integer>>> distinctDegreeFactorizatonCache;
	private Map<Polynomial<FFE>, Integer> distinctDegreeFactorizationMaxDegree;
	private Map<Polynomial<FFE>, Map<Integer, List<Polynomial<FFE>>>> irreducibleFactorizatonCache;
	private FFE primitiveRoot;

	public static class FFE extends AbstractElement<FFE> implements AlgebraicExtensionElement<PFE, FFE> {
		private UnivariatePolynomial<PFE> e;

		private FFE(FiniteField field, UnivariatePolynomial<PFE> e) {
			this.e = e;
		}

		@Override
		public int compareTo(FFE o) {
			return e.compareTo(o.e);
		}

		@Override
		public UnivariatePolynomial<PFE> asPolynomial() {
			return e;
		}

		@Override
		public String toString() {
			return this.e.toString(/*"α",*/ true);
		}

	}

	public static FiniteField getFiniteField(PrimeField prime) {
		if (finiteFields.containsKey(prime.characteristic())) {
			return finiteFields.get(prime.characteristic());
		}
		FiniteField field = new FiniteField(prime);
		finiteFields.put(prime.characteristic(), field);
		return field;
	}

	public static FiniteField getFiniteField(int q) {
		return getFiniteField(BigInteger.valueOf(q));
	}

	public static FiniteField getFiniteField(BigInteger q) {
		if (finiteFields.containsKey(q)) {
			return finiteFields.get(q);
		}
		FactorizationResult<IntE, IntE> factorization = Integers.z().uniqueFactorization(new IntE(q));
		if (factorization.primeFactors().size() != 1) {
			throw new ArithmeticException("FiniteField does not exist");
		}
		IntE prime = factorization.primeFactors().iterator().next();
		PrimeField base = PrimeField.getPrimeField(prime.getValue());
		int degree = factorization.multiplicity(prime);
		FiniteField field = getFiniteField(base, degree);
		return field;
	}

	private static <T extends Element<T>> UnivariatePolynomial<T> findIrreduciblePolynomial(Field<T> base, int degree) {
		UnivariatePolynomialRing<T> r = base.getUnivariatePolynomialRing();
		for (UnivariatePolynomial<T> result : r.monicPolynomialSet(degree)) {
			if (base.isIrreducible(result)) {
				return result;
			}
		}
		throw new ArithmeticException("No irreducible polynomials of degree " + degree + " over " + base);
	}

	public static FiniteField getFiniteField(int prime, int degree) {
		return getFiniteField(PrimeField.getPrimeField(prime), degree);
	}

	public static FiniteField getFiniteField(BigInteger prime, int degree) {
		return getFiniteField(PrimeField.getPrimeField(prime), degree);
	}

	public static FiniteField getFiniteField(PrimeField prime, int degree) {
		BigInteger q = prime.characteristic().pow(degree);
		if (finiteFields.containsKey(q)) {
			return finiteFields.get(q);
		}
		if (degree == 1) {
			return getFiniteField(prime);
		}
		FiniteField field = getFiniteField(findIrreduciblePolynomial(prime, degree), prime);
		finiteFields.put(prime.characteristic().pow(degree), field);
		return field;
	}

	public static FiniteField getFiniteField(FiniteField base, int degree) {
		return getFiniteField(findIrreduciblePolynomial(base, degree), base);
	}

	public static FiniteField getFiniteField(UnivariatePolynomial<PFE> minimalPolynomial, PrimeField prime) {
		BigInteger q = prime.characteristic().pow(minimalPolynomial.degree());
		if (finiteFields.containsKey(q)) {
			FiniteField field = finiteFields.get(q);
			if (field.minimalPolynomial().equals(minimalPolynomial)) {
				return field;
			} else {
				return new FiniteField(minimalPolynomial, prime);
			}
		}
		FiniteField field = new FiniteField(minimalPolynomial, prime);
		finiteFields.put(q, field);
		return field;
	}

	public static FiniteField getFiniteField(UnivariatePolynomial<FFE> minimalPolynomial, FiniteField base) {
		return (FiniteField) base.getExtension(minimalPolynomial).extension();
	}

	private FiniteField(PrimeField base) {
		super(base);
		init();
	}

	private FiniteField(UnivariatePolynomial<PFE> minimalpolynomial, PrimeField base) {
		super(minimalpolynomial, base, "x");
		init();
	}

	private void init() {
		this.distinctDegreeFactorizatonCache = new TreeMap<>();
		this.distinctDegreeFactorizationMaxDegree = new TreeMap<>();
		this.irreducibleFactorizatonCache = new TreeMap<>();
	}

	@Override
	public FiniteField makeExtension(UnivariatePolynomial<PFE> minimalPolynomial) {
		return getFiniteField(minimalPolynomial, PrimeField.getPrimeField(characteristic()));
	}

	public FieldEmbedding<PFE, FFE, FiniteField> getEmbeddedExtension(int degree) {
		return getEmbeddedExtension(findIrreduciblePolynomial(this, degree));
	}

	@Override
	public PrimeField getBaseField() {
		return PrimeField.getPrimeField(characteristic());
	}

	@Override
	public FFE normOver(FFE t, FieldEmbedding<PFE, FFE, FiniteField> base) {
		if (!base.getField().equals(this)) {
			throw new ArithmeticException("Not embedded here!");
		}
		BigInteger q1 = getNumberOfElements().subtract(BigInteger.ONE);
		BigInteger p1 = base.getEmbeddedField().getNumberOfElements().subtract(BigInteger.ONE);
		FFE normEmbedded = power(t, q1.divide(p1));
		UnivariatePolynomial<FFE> normPolynomial = base.asPolynomial(normEmbedded);
		if (normPolynomial.degree() > 0) {
			throw new ArithmeticException("norm not computed correctly");
		}
		return normPolynomial.univariateCoefficient(0);
	}

	@Override
	public PFE norm(FFE t) {
		BigInteger q1 = getNumberOfElements().subtract(BigInteger.ONE);
		BigInteger p1 = characteristic().subtract(BigInteger.ONE);
		FFE normEmbedded = power(t, q1.divide(p1));
		UnivariatePolynomial<PFE> normPolynomial = asPolynomial(normEmbedded);
		if (normPolynomial.degree() > 0) {
			throw new ArithmeticException("norm not computed correctly");
		}
		return normPolynomial.univariateCoefficient(0);
	}

	@Override
	public boolean hasCharacteristicRoot(FFE t, int power) {
		return true;
	}

	@Override
	public FFE characteristicRoot(FFE t, int power) {
		return power(t, characteristic().pow(degree() - (power % degree())));
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<FFE> t) {
		BigInteger q = getNumberOfElements();
		UnivariatePolynomialRing<FFE> ring = getUnivariatePolynomialRing();
		GenericAlgebraicRingExtension<FFE> cr = new GenericAlgebraicRingExtension<>(t, this);
		GenericAlgebraicExtensionElement<FFE> x = cr.fromPolynomial(ring.getVar());
		GenericAlgebraicExtensionElement<FFE> xqi = x;
		int degree = t.degree();
		for (int i = 1; i < degree; i++) {
			xqi = cr.power(xqi, q);
			GenericAlgebraicExtensionElement<FFE> xqimx = cr.subtract(xqi, x);
			if (ring.gcd(xqimx.asPolynomial(), t).degree() != 0) {
				return false;
			}
		}
		return true;
	}

	@Override
	public FactorizationResult<Polynomial<FFE>, FFE> factorization(UnivariatePolynomial<FFE> t) {
		UnivariatePolynomialRing<FFE> ring = this.getUnivariatePolynomialRing();
		SortedMap<Polynomial<FFE>, Integer> result = new TreeMap<>();
		FFE unit = t.leadingCoefficient();
		t = ring.normalize(t);
		FactorizationResult<Polynomial<FFE>, FFE> squareFreeFactors = ring.squareFreeFactorization(t);
		for (Polynomial<FFE> sff : squareFreeFactors.primeFactors()) {
			int power = squareFreeFactors.multiplicity(sff);
			sff = ring.normalize(sff);
			for (Polynomial<FFE> factor : factorizeSquareFree(sff)) {
				result.put(factor, power);
			}
		}
		return new FactorizationResult<>(unit, result);
	}

	private List<Polynomial<FFE>> factorizeSquareFree(Polynomial<FFE> t) {
		List<Polynomial<FFE>> result = new ArrayList<>();
		List<Pair<Polynomial<FFE>, Integer>> distinctDegreeFactors = this.distinctDegreeFactorization(t, -1);
		for (Pair<Polynomial<FFE>, Integer> ddf : distinctDegreeFactors) {
			result.addAll(this.irreducibleFactorization(ddf.getFirst(), ddf.getSecond()));
		}
		return result;
	}

	/**
	 * Factorizes a square free polynomial into factors that of a distinct degree.
	 * 
	 * @param t         univariate, normalized, square free polynomial over a finite
	 *                  field.
	 * @param maxDegree maximum degree of results needed. Set to -1 for all factors.
	 * @return list of pairs of polynomials and integers, such that the polynomials
	 *         are factors of t, that each consists of irreducible factors of the
	 *         the given degree.
	 */
	public List<Pair<Polynomial<FFE>, Integer>> distinctDegreeFactorization(Polynomial<FFE> t, int maxDegree) {
		if (this.distinctDegreeFactorizatonCache.containsKey(t)) {
			int previousMaxDegree = this.distinctDegreeFactorizationMaxDegree.get(t);
			if (previousMaxDegree == -1 || previousMaxDegree >= maxDegree) {
				return this.distinctDegreeFactorizatonCache.get(t);
			}
		}
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Not univariate.");
		}
		if (!t.leadingCoefficient().equals(this.one())) {
			throw new RuntimeException("Not normalized.");
		}
		UnivariatePolynomialRing<FFE> ring = this.getUnivariatePolynomialRing();
		List<Pair<Polynomial<FFE>, Integer>> result = new ArrayList<>();
		UnivariatePolynomial<FFE> work = ring.toUnivariate(t);
		GenericAlgebraicRingExtension<FFE> cr = new GenericAlgebraicRingExtension<>(work, this);
		GenericAlgebraicExtensionElement<FFE> x = cr.fromPolynomial(ring.getVar());
		GenericAlgebraicExtensionElement<FFE> xqi = x;
		for (int degree = 1; work.degree() >= 2 * degree; degree++) {
			if (maxDegree >= 0 && degree > maxDegree) {
				break;
			}
			xqi = cr.power(xqi, this.getNumberOfElements());
			GenericAlgebraicExtensionElement<FFE> xqimx = cr.subtract(xqi, x);
			Polynomial<FFE> gcd = ring.normalize(ring.gcd(work, xqimx.asPolynomial()));
			if (!gcd.equals(ring.one())) {
				result.add(new Pair<>(gcd, degree));
				work = ring.toUnivariate(ring.divideChecked(work, gcd));
				cr = new GenericAlgebraicRingExtension<FFE>(work, this);
				x = cr.fromPolynomial(ring.getVar());
				xqi = x;
				for (int j = 0; j < degree; j++) {
					xqi = cr.power(xqi, this.getNumberOfElements());
				}
			}
		}
		if (!work.equals(ring.one())) {
			result.add(new Pair<>(work, work.degree()));
		}
		this.distinctDegreeFactorizatonCache.put(t, result);
		this.distinctDegreeFactorizationMaxDegree.put(t, maxDegree);
		return result;
	}

	private List<Polynomial<FFE>> irreducibleFactorization(Polynomial<FFE> t, int degree) {
		if (!this.irreducibleFactorizatonCache.containsKey(t)) {
			this.irreducibleFactorizatonCache.put(t, new TreeMap<>());
		}
		Map<Integer, List<Polynomial<FFE>>> cache = this.irreducibleFactorizatonCache.get(t);
		if (cache.containsKey(degree)) {
			return cache.get(degree);
		}
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Not univariate.");
		}
		if (!t.leadingCoefficient().equals(this.one())) {
			throw new RuntimeException("Not normalized.");
		}
		if (t.degree() % degree != 0) {
			throw new RuntimeException("Inputs impossible!");
		}
		UnivariatePolynomialRing<FFE> ring = this.getUnivariatePolynomialRing();
		int numberOfFactors = t.degree() / degree;
		BigInteger one = BigInteger.ONE;
		BigInteger exponent = this.getNumberOfElements().pow(degree).subtract(one).shiftRight(1);
		GenericAlgebraicRingExtension<FFE> cr = new GenericAlgebraicRingExtension<>(ring.toUnivariate(t), this);
		List<Polynomial<FFE>> factors = new ArrayList<>();
		factors.add(t);
		while (factors.size() < numberOfFactors) {
			// Generate random polynomial with degree smaller than t.
			UnivariatePolynomial<FFE> h;
			do {
				h = ring.zero();
				for (int i = 0; i < t.degree(); i++) {
					h = ring.add(h, ring.multiply(this.getRandomElement(), ring.getVarPower(1, i)));
				}
			} while (h.equals(ring.zero()) || h.equals(ring.one()) || h.equals(ring.negative(ring.one())));
			Polynomial<FFE> g;
			if (characteristic().equals(BigInteger.TWO)) {
				g = h;
				GenericAlgebraicExtensionElement<FFE> hEmbedded = cr.fromPolynomial(h);
				for (int i = 0; i < this.degree() * degree - 1; i++) {
					hEmbedded = cr.multiply(hEmbedded, hEmbedded);
					g = ring.add(g, hEmbedded.asPolynomial());
				}
			} else {
				g = cr.subtract(cr.power(cr.fromPolynomial(h), exponent), cr.one()).asPolynomial();
			}
			List<Polynomial<FFE>> newFactors = new ArrayList<Polynomial<FFE>>();
			for (Polynomial<FFE> factor : factors) {
				if (factor.degree() == degree) {
					newFactors.add(factor);
					continue;
				}
				Polynomial<FFE> gcd = ring.normalize(ring.gcd(factor, g));
				if (!gcd.equals(ring.one()) && !gcd.equals(factor)) {
					newFactors.add(gcd);
					newFactors.add(ring.divideChecked(factor, gcd));
				} else {
					Polynomial<FFE> gcd2 = ring.normalize(ring.gcd(factor, ring.add(g, ring.one())));
					if (!gcd2.equals(ring.one()) && !gcd2.equals(factor)) {
						newFactors.add(gcd2);
						newFactors.add(ring.divideChecked(factor, gcd2));
					} else {
						newFactors.add(factor);
					}
				}
			}
			factors = newFactors;
		}
		cache.put(degree, factors);
		return factors;
	}

	@Override
	public boolean hasRoots(Polynomial<FFE> t) {
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Not implemented.");
		}
		UnivariatePolynomialRing<FFE> ring = this.getUnivariatePolynomialRing();
		GenericAlgebraicRingExtension<FFE> cr = new GenericAlgebraicRingExtension<>(ring.toUnivariate(t), this);
		GenericAlgebraicExtensionElement<FFE> x = cr.fromPolynomial(ring.getVar());
		GenericAlgebraicExtensionElement<FFE> xq = cr.power(x, this.getNumberOfElements());
		Polynomial<FFE> xqmx = cr.subtract(xq, x).asPolynomial();
		Polynomial<FFE> gcd = ring.gcd(t, xqmx);
		return gcd.degree() > 0;
	}

	@Override
	public Map<FFE, Integer> roots(Polynomial<FFE> t) {
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Multivariate root finding not implemented. Use ideal logic.");
		}
		UnivariatePolynomialRing<FFE> ring = getUnivariatePolynomialRing();
		Map<FFE, Integer> result = new TreeMap<>();
		t = ring.normalize(t);
		FactorizationResult<Polynomial<FFE>, FFE> squareFreeFactors = ring.squareFreeFactorization(t);
		for (Polynomial<FFE> sff : squareFreeFactors.primeFactors()) {
			sff = ring.normalize(sff);
			List<Pair<Polynomial<FFE>, Integer>> distinctDegreeFactors = this.distinctDegreeFactorization(sff, 1);
			for (Pair<Polynomial<FFE>, Integer> ddf : distinctDegreeFactors) {
				if (ddf.getSecond() != 1) {
					continue;
				}
				List<Polynomial<FFE>> linear = this.irreducibleFactorization(ddf.getFirst(), ddf.getSecond());
				for (Polynomial<FFE> linearFactor : linear) {
					result.put(this.negative(linearFactor.coefficient(ring.getMonomial(new int[] { 0 }))),
							squareFreeFactors.multiplicity(sff));
				}
			}
		}
		return result;
	}

	@Override
	public FFE primitiveRoot() {
		if (primitiveRoot != null) {
			return primitiveRoot;
		}
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(new IntE(getNumberOfUnits()));
		List<BigInteger> tests = new ArrayList<>();
		for (IntE prime : factors.primeFactors()) {
			tests.add(getNumberOfUnits().divide(prime.getValue()));
		}
		candidateLoop: for (FFE candidate : this) {
			candidate = add(alpha(), candidate);
			if (candidate.equals(zero())) {
				continue;
			}
			for (BigInteger test : tests) {
				if (power(candidate, test).equals(one())) {
					continue candidateLoop;
				}
			}
			this.primitiveRoot = candidate;
			break;
		}
		return primitiveRoot;
	}

	public IntE discreteLogarithm(FFE base, FFE power) {
		Integers z = Integers.z();
		IntE numberOfUnits = z.getInteger(getNumberOfUnits());
		FFE primitiveRoot = primitiveRoot();
		if (!base.equals(primitiveRoot)) {
			ModuloIntegerRing mod = new ModuloIntegerRing(numberOfUnits.getValue());
			IntE powerLog = discreteLogarithm(primitiveRoot, power);
			IntE baseLog = discreteLogarithm(primitiveRoot, base);
			IntE gcd = z.gcd(powerLog, baseLog);
			powerLog = z.divideChecked(powerLog, gcd);
			baseLog = z.divideChecked(baseLog, gcd);
			return z.multiply(powerLog, mod.lift(mod.inverse(mod.reduce(baseLog))));
		}
		NumberFieldIntegers nf = NumberField.getNumberFieldFromIntegerPolynomial(
				z.centeredLiftUnivariatePolynomial(minimalPolynomial(), characteristic())).maximalOrder();
		IntE smoothnessBound = z.getInteger(getNumberOfElements().bitLength() + 1);
		List<NumberFieldIdeal> primeIdeals = new ArrayList<>();
		Set<IntE> smoothnessBase = new TreeSet<>();
		for (IntE prime : z.setOfPrimes()) {
			if (prime.compareTo(smoothnessBound) > 0) {
				break;
			}
			if (z.isDivisible(nf.discriminant(), prime)) {
				continue;
			}
			smoothnessBase.add(prime);
			primeIdeals.addAll(nf.idealsOver(prime));
		}
		NumberFieldIdeal ideal = nf.idealsOver(characteristic()).get(0);
		List<List<IntE>> relations = new ArrayList<>();
		List<IntE> exponents = new ArrayList<>();
		Matrix<IntE> relationsMatrix;
		MatrixModule<IntE> relationsMatrixModule;
		while (true) {
			IntE exponent = z.getRandomElement(numberOfUnits);
			FFE testPower = power(primitiveRoot, exponent);
			NFE lifted = nf.centeredLift(testPower, ideal);
			Optional<NumberFieldIdeal> principal = nf.getIdealIfSmoothOver(Collections.singletonList(lifted),
					smoothnessBase);
			if (principal.isEmpty()) {
				continue;
			}
			FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = nf.idealFactorization(principal.get());
			exponents.add(exponent);
			List<IntE> relation = new ArrayList<>();
			relation.addAll(nf.asUnitGeneratorVector(nf.projectToUnit(lifted)).asList());
			for (NumberFieldIdeal primeIdeal : primeIdeals) {
				relation.add(z.getInteger(factorization.multiplicity(primeIdeal)));
			}
			relations.add(relation);
			relationsMatrix = new Matrix<>(relations);
			relationsMatrixModule = relationsMatrix.getModule(z);
			if (relationsMatrixModule.rank(relationsMatrix) == relationsMatrix.columns()) {
				break;
			}
		}
		Vector<IntE> logarithms = relationsMatrixModule.solve(relationsMatrix, new Vector<>(exponents));
		while (true) {
			IntE exponent = z.getRandomElement(z.getInteger(getNumberOfUnits()));
			FFE testPower = multiply(power, power(primitiveRoot, z.negative(exponent)));
			NFE lifted = nf.centeredLift(testPower, ideal);
			Optional<NumberFieldIdeal> principal = nf.getIdealIfSmoothOver(Collections.singletonList(lifted),
					smoothnessBase);
			if (principal.isEmpty()) {
				continue;
			}
			FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = nf.idealFactorization(principal.get());
			IntE result = exponent;
			Vector<IntE> asUnitVector = nf.asUnitGeneratorVector(nf.projectToUnit(lifted));
			for (int i = 0; i < asUnitVector.dimension(); i++) {
				result = z.add(z.multiply(logarithms.get(i + 1), asUnitVector.get(i + 1)), result);
				result = z.remainder(result, numberOfUnits);
			}
			for (int i = 0; i < primeIdeals.size(); i++) {
				NumberFieldIdeal primeIdeal = primeIdeals.get(i);
				IntE multiplicity = z.getInteger(factorization.multiplicity(primeIdeal));
				result = z.add(z.multiply(logarithms.get(asUnitVector.dimension() + i + 1), multiplicity), result);
				result = z.remainder(result, numberOfUnits);
			}
			return result;
		}
//		System.out.println(logarithms);
//
//		Pair<UnivariatePolynomial<IntE>, UnivariatePolynomial<IntE>> polynomialChoice = choosePolynomial(degree() + 2);
//		return null;
	}

//	private Pair<UnivariatePolynomial<IntE>, UnivariatePolynomial<IntE>> choosePolynomial(int d) {
//		Integers z = Integers.z();
//		PrimeField base = PrimeField.getPrimeField(characteristic());
//		UnivariatePolynomialRing<IntE> integerPolynomials = z.getUnivariatePolynomialRing();
//		UnivariatePolynomialRing<PFE> polynomials = base.getUnivariatePolynomialRing();
//		UnivariatePolynomial<PFE> psi = minimalPolynomial();
//		UnivariatePolynomial<IntE> quadratic = null;
//		PFE root = null;
//		Optional<Fraction> reconstruction = Optional.empty();
//		for (UnivariatePolynomial<IntE> q : integerPolynomials.monicPolynomialSet(2)) {
//			Map<PFE, Integer> roots = base.roots(z.reduceUnivariatePolynomial(q, characteristic()));
//			if (roots.size() != 2) {
//				continue;
//			}
//			Iterator<PFE> it = roots.keySet().iterator();
//			root = it.next();
//			if (root.equals(base.zero())) {
//				root = it.next();
//			}
//			reconstruction = base.rationalReconstruction(root);
//			if (reconstruction.isEmpty()) {
//				continue;
//			}
//			if (!z.factorization(q).isIrreducible()) {
//				continue;
//			}
//			quadratic = q;
//			break;
//		}
//		UnivariatePolynomial<IntE> f = null;
//		UnivariatePolynomial<IntE> g = null;
//		UnivariatePolynomial<IntE> g0 = null;
//		UnivariatePolynomial<IntE> g1 = null;
//		PolynomialRing<IntE> int2 = AbstractPolynomialRing.getPolynomialRing(z, 2, Monomial.GREVLEX);
//		for (UnivariatePolynomial<IntE> g1c : integerPolynomials.polynomialSet(degree() - 1)) {
//			g1 = g1c;
//			if (g1.degree() < 1) {
//				continue;
//			}
//			g0 = z.centeredLiftUnivariatePolynomial(
//					polynomials.subtract(psi,
//							polynomials.multiply(root, z.reduceUnivariatePolynomial(g1, characteristic()))),
//					characteristic());
//			g = integerPolynomials.add(integerPolynomials.scalarMultiply(reconstruction.get().getDenominator(), g0),
//					integerPolynomials.scalarMultiply(reconstruction.get().getNumerator(), g1));
//			if (!z.factorization(g).isIrreducible()) {
//				continue;
//			}
//			f = integerPolynomials.toUnivariate(int2.resultant(int2.getEmbedding(quadratic, new int[] { 1 }),
//					int2.add(int2.getEmbedding(g0, new int[] { 0 }),
//							int2.multiply(int2.getEmbedding(g1, new int[] { 0 }), int2.getVar(2))),
//					2));
//			if (!z.factorization(f).isIrreducible()) {
//				continue;
//			}
//			break;
//		}
////		UnivariatePolynomial<IntE> f = null;
////		for (UnivariatePolynomial<IntE> q : integerPolynomials.monicPolynomialSet(d + 1 - degree())) {
////			f = z.centeredLiftUnivariatePolynomial(
////					polynomials.multiply(z.reduceUnivariatePolynomial(q, characteristic()), psi),
////					base.characteristic());
////			if (z.factorization(f).isIrreducible()) {
////				break;
////			}
////		}
////		List<Vector<IntE>> basis = new ArrayList<>();
////		for (int i = 0; i < degree(); i++) {
////			List<IntE> c = new ArrayList<>();
////			for (int j = 0; j < d + 1; j++) {
////				c.add(z.zero());
////			}
////			c.set(i, z.getInteger(characteristic()));
////			basis.add(new Vector<>(c));
////		}
////		for (int i = degree(); i < d + 1; i++) {
////			List<IntE> c = new ArrayList<>();
////			for (int j = 0; j < d + 1; j++) {
////				int index = degree() - i + j;
////				if (index >= 0 && index <= degree()) {
////					c.add(z.centeredLift(psi.univariateCoefficient(index), characteristic()));
////				} else {
////					c.add(z.zero());
////				}
////			}
////			basis.add(new Vector<>(c));
////		}
////		basis = new FiniteRationalVectorSpace(d + 1).integerLatticeReduction(basis);
////		List<IntE> c = new ArrayList<>();
////		for (Vector<IntE> b : basis) {
////			c.add(b.get(1));
////		}
////		UnivariatePolynomial<IntE> g = integerPolynomials.getPolynomial(c);
//		if (!f.leadingCoefficient().equals(z.one()) || !z.factorization(f).isIrreducible()
//				|| !z.factorization(g).isIrreducible()
//				|| !polynomials.gcd(z.reduceUnivariatePolynomial(f, characteristic()),
//						z.reduceUnivariatePolynomial(g, characteristic())).equals(psi)) {
//			throw new ArithmeticException("Polynomial choice wrong");
//		}
//		return new Pair<>(f, g);
//	}

	@Override
	public GaloisGroup<PFE, FFE, FiniteField> galoisGroup() {
		return new GaloisGroup<>(this, cyclicGenerator());
	}

	@Override
	public FFE fromSmallDegreePolynomial(UnivariatePolynomial<PFE> t) {
		return new FFE(this, t);
	}

	@Override
	public List<FFE> conjugates(FFE s) {
		List<FFE> conjugates = new ArrayList<>();
		for (int i = 0; i < degree(); i++) {
			conjugates.add(power(s, getBaseField().getNumberOfElements().pow(i)));
		}
		return conjugates;
	}

	@Override
	public boolean isCyclic() {
		return true;
	}

	@Override
	public FieldAutomorphism<PFE, FFE, FiniteField> cyclicGenerator() {
		int[] images = new int[degree()];
		for (int i = 0; i < degree(); i++) {
			images[i] = (i + 1) % degree();
		}
		return new FieldAutomorphism<>(this, images);
	}

	public int smallestSubfieldContaining(FFE t) {
		FFE power = t;
		int result = 0;
		while (true) {
			result++;
			power = power(power, characteristic());
			if (power.equals(t)) {
				return result;
			}
		}
	}

	@Override
	protected FiniteField asExtensionType() {
		return this;
	}
}
