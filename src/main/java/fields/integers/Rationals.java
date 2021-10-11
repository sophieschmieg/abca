package fields.integers;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.helper.FieldOfFractions;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.polynomials.Monomial;
import fields.vectors.FreeModule;
import fields.vectors.Vector;
import util.MiscAlgorithms;

public class Rationals extends AbstractField<Fraction> {
	private static Integers z = Integers.z();
	private static FieldOfFractions<IntE> q = new FieldOfFractions<>(z);
	private static Rationals rationals = new Rationals();
	private Fraction zero;
	private Fraction one;
	private Map<BigInteger, LocalizedFractions> localized;
	private ValueFractions infValue;

	public static class Fraction extends AbstractElement<Fraction> {
		private FieldOfFractions.Fraction<IntE> fraction;

		public Fraction(FieldOfFractions.Fraction<IntE> fraction) {
			this.fraction = fraction;
		}

		public FieldOfFractions.Fraction<IntE> fraction() {
			return fraction;
		}

		public void canonicalize() {
			fraction.canonicalize();
		}

		public IntE getNumerator() {
			canonicalize();
			return fraction.getNumerator();
		}

		public IntE getDenominator() {
			canonicalize();
			return fraction.getDenominator();
		}

		@Override
		public int compareTo(Fraction o) {
			return fraction.compareTo(o.fraction);
		}

		@Override
		public String toString() {
			return fraction.toString();
		}

		public int numeratorIntValueExact() {
			canonicalize();
			return getNumerator().getValue().intValueExact();
		}

		public int denominatorIntValueExact() {
			canonicalize();
			return getDenominator().getValue().intValueExact();
		}

		public IntE round() {
			BigInteger reducedNumerator = getNumerator().getValue().mod(getDenominator().getValue());
			if (reducedNumerator.shiftLeft(1).compareTo(getDenominator().getValue()) > 0) {
				return roundUp();
			}
			return roundDown();
		}

		public IntE roundUp() {
			canonicalize();
			BigInteger mod = getNumerator().getValue().mod(getDenominator().getValue());
			BigInteger roundUpNumerator = mod.equals(BigInteger.ZERO) ? getNumerator().getValue()
					: getNumerator().getValue().add(getDenominator().getValue()).subtract(mod);
			return new IntE(roundUpNumerator.divide(getDenominator().getValue()));
		}

		public IntE roundDown() {
			canonicalize();
			BigInteger mod = getNumerator().getValue().mod(getDenominator().getValue());
			BigInteger roundDownNumerator = getNumerator().getValue().subtract(mod);
			return new IntE(roundDownNumerator.divide(getDenominator().getValue()));
		}

		public IntE asInteger() {
			canonicalize();
			if (!getDenominator().equals(z.one())) {
				throw new ArithmeticException("Not an integer!");
			}
			return getNumerator();
		}
	}

	private Rationals() {
		zero = new Fraction(q.getEmbedding(z.zero()));
		one = new Fraction(q.getEmbedding(z.one()));
		localized = new TreeMap<>();
	}

	public static Rationals q() {
		return rationals;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public Extension<Fraction, Fraction, NFE, NumberField> getExtension(
			UnivariatePolynomial<Fraction> minimalPolynomial) {
		NumberField extension = new NumberField(minimalPolynomial);
		return new Extension<>(extension, this, extension.getEmbeddingMap(), extension.asVectorMap());
	}

	@Override
	public FactorizationResult<Polynomial<Fraction>, Fraction> factorization(UnivariatePolynomial<Fraction> t) {
		UnivariatePolynomialRing<Fraction> ring = this.getUnivariatePolynomialRing();
		SortedMap<Polynomial<Fraction>, Integer> result = new TreeMap<>();
		Fraction unit = t.leadingCoefficient();
		t = ring.normalize(t);
		IntE lcm = t.leadingCoefficient().getDenominator();
		for (int i = 0; i < t.degree(); i++) {
			lcm = z.lcm(lcm, t.univariateCoefficient(i).getDenominator());
		}
		t = ring.multiply(getInteger(lcm), t);
		UnivariatePolynomialRing<IntE> zRing = z.getUnivariatePolynomialRing();
		UnivariatePolynomial<IntE> asIntegerPolynomial = zRing.getEmbedding(t, new MathMap<>() {
			@Override
			public IntE evaluate(Fraction t) {
				return t.getNumerator();
			}
		});
		FactorizationResult<Polynomial<IntE>, IntE> factors = z.factorization(asIntegerPolynomial);
		for (Polynomial<IntE> factor : factors.primeFactors()) {
			if (factor.degree() <= 0) {
				continue;
			}
			result.put(ring.normalize(ring.getEmbedding(factor, new MathMap<>() {
				@Override
				public Fraction evaluate(IntE t) {
					return getInteger(t);
				}
			})), factors.multiplicity(factor));
		}
		return new FactorizationResult<>(unit, result);
	}

	@Override
	public Fraction zero() {
		return zero;
	}

	@Override
	public Fraction one() {
		return one;
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public Fraction add(Fraction t1, Fraction t2) {
		return new Fraction(q.add(t1.fraction, t2.fraction));
	}

	@Override
	public Fraction negative(Fraction t) {
		return new Fraction(q.negative(t.fraction()));
	}

	@Override
	public Fraction multiply(Fraction t1, Fraction t2) {
		return new Fraction(q.multiply(t1.fraction, t2.fraction));
	}

	@Override
	public Fraction inverse(Fraction t) {
		return new Fraction(q.inverse(t.fraction()));
	}

	public LocalizedFractions withValuation(BigInteger prime) {
		if (!localized.containsKey(prime)) {
			localized.put(prime, new LocalizedFractions(prime));
		}
		return localized.get(prime);
	}

	public ValueFractions withInfValue() {
		if (infValue == null) {
			infValue = new ValueFractions();
		}
		return infValue;
	}

	@Override
	public Fraction getRandomElement() {
		return /* getInteger(z.getRandomElement(new IntE(20))); */new Fraction(q.getRandomElement());
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
	public Iterator<Fraction> iterator() {
		return new Iterator<>() {
			Iterator<FieldOfFractions.Fraction<IntE>> it = q.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Fraction next() {
				return new Fraction(it.next());
			}
		};
	}

	public Fraction getEmbedding(IntE t) {
		return new Fraction(q.getEmbedding(t));
	}

	public MathMap<IntE, Fraction> getEmbeddingMap() {
		return new MathMap<>() {

			@Override
			public Fraction evaluate(IntE t) {
				return getInteger(t);
			}
		};
	}

	public Fraction getFraction(IntE numerator, IntE denominator) {
		return new Fraction(q.getFraction(numerator, denominator));
	}

	public Fraction getFraction(int numerator, int denominator) {
		return getFraction(new IntE(numerator), new IntE(denominator));
	}

	public Fraction getFraction(BigInteger numerator, BigInteger denominator) {
		return getFraction(new IntE(numerator), new IntE(denominator));
	}

	public Fraction getEmbedding(BigInteger t) {
		return getEmbedding(z.getInteger(t));
	}

	public Fraction getInteger(BigInteger t) {
		return getEmbedding(t);
	}

	public boolean isInteger(Fraction t) {
		return t.fraction().getDenominator().equals(z.one());
	}

	@Override
	public String toString() {
		return "Q";
	}

	public Iterator<IntE> continuedFraction(Fraction t) {
		return MiscAlgorithms.continuedFraction(this, t, new MathMap<>() {

			@Override
			public IntE evaluate(Fraction t) {
				return t.roundDown();
			}
		});
	}

	public Iterator<Fraction> continuedFractionApproximation(Fraction t) {
		return MiscAlgorithms.continuedFractionApproximation(this, t, new MathMap<>() {

			@Override
			public IntE evaluate(Fraction t) {
				return t.roundDown();
			}
		});
	}

	public Optional<Fraction> rationalReconstruction(IntE t, IntE m, BigInteger numeratorBound, BigInteger denominatorBound) {
		Integers z = Integers.z();
		FreeModule<IntE> vectors = new FreeModule<>(z, 2);
		Vector<IntE> v = new Vector<>(m, z.zero());
		Vector<IntE> w = new Vector<>(t, z.one());
		while (w.get(1).getValue().compareTo(numeratorBound) > 0) {
			IntE q = getFraction(v.get(1), w.get(1)).roundDown();
			Vector<IntE> next = vectors.subtract(v, vectors.scalarMultiply(q, w));
			v = w;
			w = next;
		}
		if (w.get(2).compareTo(z.zero()) < 0) {
			w = vectors.negative(w);
		}
		if (w.get(2).getValue().compareTo(denominatorBound) > 0 || !z.gcd(w.get(1), w.get(2)).equals(z.one())) {
			return Optional.empty();
		}
		return Optional.of(getFraction(w.get(1), w.get(2)));
	}

	public Optional<Fraction> rationalReconstruction(IntE t, IntE m) {
		if (m.equals(new IntE(2))) {
			return rationalReconstruction(t, m, BigInteger.ONE, BigInteger.ONE);
		}
		BigInteger sqrt = m.getValue().shiftRight(1).sqrt();
		return rationalReconstruction(t, m, sqrt.add(BigInteger.ONE), sqrt);
	}
	
	public Optional<Polynomial<Fraction>> reconstructPolynomial(PolynomialRing<Fraction> polynomialRing, List<Polynomial<PFE>> perPrime, List<PrimeField> fields) {
		List<Ideal<IntE>> asIdealList = new ArrayList<>();
		for (PrimeField fp : fields) {
			asIdealList.add(z.getIdeal(Collections.singletonList(new IntE(fp.getNumberOfElements()))));
		}
		ChineseRemainderPreparation<IntE> preparation = z.prepareChineseRemainderTheorem(asIdealList);
		IntE modulus = preparation.getProduct().generators().get(0);
		Set<Monomial> monomials = new TreeSet<>();
		for (Polynomial<PFE> polynomial : perPrime) {
			monomials.addAll(polynomial.monomials());
		}
		Map<Monomial, Fraction> coefficients = new TreeMap<>();
		for (Monomial m : monomials) {
			List<IntE> lifts = new ArrayList<>();
			for (int i = 0; i < perPrime.size(); i++) {
				PrimeField fp = fields.get(i);
				lifts.add(fp.liftToInteger(perPrime.get(i).coefficient(m)));
			}
			IntE crt = z.chineseRemainderTheorem(lifts, preparation);
			Optional<Fraction> reconstruction = rationalReconstruction(crt, modulus);
			if (reconstruction.isEmpty()) {
				return Optional.empty();
			}
			coefficients.put(m, reconstruction.get());
		}
		return Optional.of(polynomialRing.getPolynomial(coefficients));
	}

}
