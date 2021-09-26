package fields.finitefields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField.PFE;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import util.MiscAlgorithms;

public class PrimeField extends AbstractField<PFE> {
	private static Map<BigInteger, PrimeField> primeFields = new TreeMap<>();

	private BigInteger p;
	private boolean verySmallPrime;
	private int smallPrime;
	private PFE[] elements;
	private PFE[][] divisiontable;
	private PFE[] sqrt;
	private PFE primitiveRoot;
	private PFE quadraticNonResidue;

	public static PrimeField getPrimeField(int p) {
		return getPrimeField(BigInteger.valueOf(p));
	}

	public static PrimeField getPrimeField(IntE p) {
		return getPrimeField(p.getValue());
	}

	public static PrimeField getPrimeField(BigInteger p) {
		if (primeFields.containsKey(p)) {
			return primeFields.get(p);
		}
		PrimeField field = new PrimeField(p);
		primeFields.put(p, field);
		return field;
	}

	private PrimeField(BigInteger p) {
		this.p = p;
		if (!p.isProbablePrime(100)) {
			throw new RuntimeException("Not a prime number");
		}
		if (p.compareTo(BigInteger.ZERO) < 0) {
			throw new RuntimeException("Not a positive prime number");
		}
		if (p.compareTo(BigInteger.valueOf(257)) > 0) {
			return;
		}
		this.smallPrime = p.intValueExact();
		this.verySmallPrime = true;
		this.elements = new PFE[smallPrime];
		this.divisiontable = new PFE[smallPrime][smallPrime];
		this.sqrt = new PFE[smallPrime];
		for (int i = 0; i < smallPrime; i++) {
			this.elements[i] = new PFE(i);
			if (this.sqrt[(i * i) % smallPrime] == null) {
				this.sqrt[(i * i) % smallPrime] = this.elements[i];
			}
		}
		/*
		 * for (int i = 0; i < smallPrime; i++) { BigInteger divisor =
		 * BigInteger.valueOf(i); for (int j = 0; j < smallPrime; j++) { if (i != 0) {
		 * this.divisiontable[j][i] =
		 * this.elements[(divisor.modInverse(p).intValueExact() * j) % smallPrime]; } }
		 * if (this.sqrt[(i * i) % smallPrime] == null) { this.sqrt[(i * i) %
		 * smallPrime] = this.elements[i]; } }
		 */
	}

	private PrimeField(int p) {
		this(BigInteger.valueOf(p));
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public PFE zero() {
		return this.getElement(0);
	}

	@Override
	public PFE one() {
		return this.getElement(1);
	}

	@Override
	public BigInteger characteristic() {
		return p;
	}

	@Override
	public PFE primitiveRoot() {
		if (primitiveRoot != null) {
			return primitiveRoot;
		}
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(new IntE(getNumberOfUnits()));
		List<BigInteger> tests = new ArrayList<>();
		for (IntE prime : factors.primeFactors()) {
			tests.add(getNumberOfUnits().divide(prime.getValue()));
		}
		candidateLoop: for (PFE candidate : this) {
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

	@Override
	public PFE add(PFE t1, PFE t2) {
		if (verySmallPrime) {
			return elements[(t1.index + t2.index) % smallPrime];
		}
		return getElement(t1.value.add(t2.value));
	}

	@Override
	public PFE negative(PFE t) {
		if (verySmallPrime) {
			return elements[t.index == 0 ? 0 : (smallPrime - t.index)];
		}
		return this.getElement(t.value.negate());
	}

	@Override
	public PFE multiply(PFE t1, PFE t2) {
		if (verySmallPrime) {
			return elements[(t1.index * t2.index) % smallPrime];
		}
		return getElement(t1.value.multiply(t2.value));
	}

	@Override
	public PFE inverse(PFE t) {
		if (verySmallPrime) {
			if (t.index == 0) {
				throw new ArithmeticException("Div by 0");
			}
			if (divisiontable[1][t.index] == null) {
				divisiontable[1][t.index] = this.elements[t.value.modInverse(p).intValueExact()];
			}
			return divisiontable[1][t.index];
		}
		return this.getElement(t.value.modInverse(p));
	}

	@Override
	public PFE divide(PFE dividend, PFE divisor) {
		if (verySmallPrime) {
			if (divisor.index == 0) {
				throw new ArithmeticException("Div by 0");
			}
			if (divisiontable[dividend.index][divisor.index] == null) {
				divisiontable[dividend.index][divisor.index] = this.elements[(divisor.value.modInverse(p)
						.intValueExact() * dividend.index) % smallPrime];
			}
			return divisiontable[dividend.index][divisor.index];
		}
		return this.getElement(dividend.value.multiply(divisor.value.modInverse(p)));
	}

	
	@Override
	public boolean hasSqrt(PFE t) {
		if (verySmallPrime) {
			return sqrt[t.index] != null;
		}
		if (this.p == BigInteger.TWO) {
			return true;
		}
		int jacobi = MiscAlgorithms.jacobiSymbol(t.value, this.p);
		return jacobi == 1 || jacobi == 0;
	}
	
	public PFE quadraticNonResidue() {
		if (p == BigInteger.TWO) {
			throw new ArithmeticException("No quadratic non residues mod 2");
		}
		if (quadraticNonResidue == null) {
			for (PFE t : this) {
				if (MiscAlgorithms.jacobiSymbol(t.value, p) != -1) {
					continue;
				}
				quadraticNonResidue = t;
				break;
			}
		}
		return quadraticNonResidue;
	}

	@Override
	public Map<PFE, Integer> sqrt(PFE t) {
		if (t.equals(this.zero())) {
			return Collections.singletonMap(this.zero(), 2);
		}
		if (verySmallPrime) {
			if (sqrt[t.index] == null) {
				return Collections.emptyMap();
			}
			if (characteristic().equals(BigInteger.TWO)) {
				return Collections.singletonMap(sqrt[t.index], 2);
			}
			Map<PFE, Integer> result = new TreeMap<>();
			result.put(sqrt[t.index], 1);
			result.put(elements[smallPrime - sqrt[t.index].index], 1);
			return result;
		}
		if (MiscAlgorithms.jacobiSymbol(t.value, this.p) == -1) {
			return Collections.emptyMap();
		}
		BigInteger zero = BigInteger.ZERO;
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		BigInteger four = BigInteger.valueOf(4);
		if (this.p.mod(four).intValue() == 3) {
			PFE sqrt = this.power(t, this.p.add(one).divide(four));
			Map<PFE, Integer> result = new TreeMap<>();
			result.put(sqrt, 1);
			result.put(this.negative(sqrt), 1);
			return result;
		}
		if (this.p.equals(two)) {
			return Collections.singletonMap(t, 2);
		}
		int s = 0;
		BigInteger q = this.p.subtract(one);
		while (q.mod(two).equals(zero)) {
			s++;
			q = q.shiftRight(1);
		}
		BigInteger z = two;
		while (z.mod(this.p).equals(zero) || MiscAlgorithms.jacobiSymbol(z, this.p) == 1) {
			z = z.add(one);
		}
		int m = s;
		PFE c = this.power(this.getElement(z), q);
		PFE n = this.power(t, q);
		PFE r = this.power(t, q.add(one).shiftRight(1)); // (q + 1) / 2
		while (!n.equals(this.one())) {
			int i = 0;
			PFE sq = n;
			while (!sq.equals(this.one())) {
				sq = this.multiply(sq, sq);
				i++;
			}
			PFE b = c;
			for (int j = 0; j < m - i - 1; j++) {
				b = this.multiply(b, b);
			}
			m = i;
			c = this.multiply(b, b);
			n = this.multiply(n, c);
			r = this.multiply(r, b);
		}
		Map<PFE, Integer> result = new TreeMap<>();
		result.put(r, 1);
		result.put(this.negative(r), 1);
		return result;
	}

	@Override
	public PFE getRandomElement() {
		return getElement(MiscAlgorithms.randomBigInteger(new Random(), p));
	}

	@Override
	public BigInteger getNumberOfElements() {
		return this.p;
	}

	@Override
	public Iterator<PFE> iterator() {
		return new Iterator<>() {
			private BigInteger i = BigInteger.ZERO;

			@Override
			public boolean hasNext() {
				return this.i.compareTo(p) < 0;
			}

			@Override
			public PFE next() {
				BigInteger tmp = this.i;
				this.i = this.i.add(BigInteger.ONE);
				return PrimeField.this.getElement(tmp);
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};

	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public String toString() {
		return "F" + p;
	}

	@Override
	public Extension<PFE, PFE, FFE, FiniteField> getExtension(UnivariatePolynomial<PFE> minimalPolynomial) {
		FiniteField extension = FiniteField.getFiniteField(minimalPolynomial, this);
		return new Extension<>(extension, this, extension.getEmbeddingMap(), extension.asVectorMap());
	}

	public PFE getElement(int value) {
		if (verySmallPrime) {
			if (value < 0) {
				return elements[smallPrime - (-value % smallPrime)];
			}
			return elements[value % smallPrime];
		}
		return new PFE(BigInteger.valueOf(value));
	}

	public PFE getElement(BigInteger value) {
		if (verySmallPrime) {
			return getElement(value.mod(p).intValueExact());
		}
		return new PFE(value);
	}

	public PFE getInteger(int value) {
		return getElement(value);
	}

	public PFE getInteger(BigInteger value) {
		return getElement(value);
	}

	public PFE reduce(IntE t) {
		return getElement(t.getValue());
	}

	public PFE reduce(Fraction t) {
		return divide(reduce(t.getNumerator()), reduce(t.getDenominator()));
	}

	public IntE lift(PFE t) {
		return Integers.z().getInteger(t.value);
	}

	public Optional<Fraction> rationalReconstruction(PFE t, BigInteger numeratorBound, BigInteger denominatorBound) {
		return Rationals.q().rationalReconstruction(lift(t), new IntE(p), numeratorBound, denominatorBound);
	}

	public Optional<Fraction> rationalReconstruction(PFE t) {
		return Rationals.q().rationalReconstruction(lift(t), new IntE(p));
	}

	@Override
	public boolean hasCharacteristicRoot(PFE t, int power) {
		return true;
	}

	@Override
	public PFE characteristicRoot(PFE t, int power) {
		return t;
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<PFE> t) {
		FiniteField base = FiniteField.getFiniteField(this);
		UnivariatePolynomialRing<FFE> baseRing = base.getUnivariatePolynomialRing();
		return base.isIrreducible(baseRing.getEmbedding(t, base.getEmbeddingMap()));
	}

	@Override
	public FactorizationResult<Polynomial<PFE>, PFE> factorization(UnivariatePolynomial<PFE> t) {
		FiniteField base = FiniteField.getFiniteField(this);
		UnivariatePolynomialRing<FFE> baseRing = base.getUnivariatePolynomialRing();
		FactorizationResult<Polynomial<FFE>, FFE> factors = base
				.factorization(baseRing.getEmbedding(t, base.getEmbeddingMap()));
		SortedMap<Polynomial<PFE>, Integer> result = new TreeMap<>();
		for (Polynomial<FFE> factor : factors.primeFactors()) {
			result.put(this.getUnivariatePolynomialRing().getEmbedding(factor, base.asBaseFieldElementMap()),
					factors.multiplicity(factor));
		}
		return new FactorizationResult<>(base.asBaseFieldElement(factors.getUnit()), result);
	}

	@Override
	public boolean hasRoots(Polynomial<PFE> t) {
		FiniteField base = FiniteField.getFiniteField(this);
		return base.hasRoots(base.getUnivariatePolynomialRing().getEmbedding(t, base.getEmbeddingMap()));
	}

	@Override
	public Map<PFE, Integer> roots(Polynomial<PFE> t) {
		FiniteField base = FiniteField.getFiniteField(this);
		Map<FFE, Integer> roots = base
				.roots(base.getUnivariatePolynomialRing().getEmbedding(t, base.getEmbeddingMap()));
		Map<PFE, Integer> result = new TreeMap<>();
		for (FFE root : roots.keySet()) {
			result.put(base.asBaseFieldElement(root), roots.get(root));
		}
		return result;
	}

	public class PFE extends AbstractElement<PFE> {
		private BigInteger value;
		private int index;

		private PFE(BigInteger value) {
			this.value = value.mod(p);
			this.index = -1;
			while (this.value.compareTo(BigInteger.ZERO) < 0)
				this.value = this.value.add(p);
		}

		private PFE(int index) {
			this.value = BigInteger.valueOf(index);
			this.index = index;
			while (this.value.compareTo(BigInteger.ZERO) < 0)
				this.value = this.value.add(p);
		}

		public BigInteger getValue() {
			return value;
		}

		public String toString() {
			return Integers.z().centeredLift(this, p).toString();
		}

		@Override
		public int compareTo(PFE o) {
			return this.value.compareTo(o.value);
		}

	}
}
