package fields.helper;

import java.io.IOException;
import java.math.BigInteger;
import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.FieldOfFractions.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import util.PeekableReader;

public class FieldOfFractions<T extends Element<T>> extends AbstractField<Fraction<T>> implements Field<Fraction<T>> {
	private Ring<T> ring;

	public FieldOfFractions(Ring<T> ring) {
		this.ring = ring;
		if (!this.ring.isIntegral())
			throw new RuntimeException("Object does not exist");
	}

	@Override
	public String toString() {
		return "Q(" + ring + ")";
	}

	@Override
	public Fraction<T> parse(PeekableReader reader) throws IOException {
		T numerator = ring.parse(reader);
		T denominator = ring.one();
		if (reader.peek() == '/') {
			reader.read();
			denominator = ring.parse(reader);
		}
		return getFraction(numerator, denominator);
	}

	@Override
	public Exactness exactness() {
		return ring.exactness();
	}

	public Ring<T> getRing() {
		return ring;
	}

	@Override
	public boolean isFinite() {
		return this.ring.isFinite();
	}

	@Override
	public Fraction<T> zero() {
		return this.getEmbedding(this.ring.zero());
	}

	@Override
	public Fraction<T> one() {
		return this.getEmbedding(this.ring.one());
	}

	@Override
	public BigInteger characteristic() {
		return this.ring.characteristic();
	}

	@Override
	public Fraction<T> add(Fraction<T> t1, Fraction<T> t2) {
		if (t1.denominator.equals(t2.denominator)) {
			return new Fraction<>(this.ring, this.ring.add(t1.numerator, t2.numerator), t1.denominator, t1.counter + t2.counter + 1);
		}
		if (t1.counter > 1000) {
			t1.canonicalize();
		}
		if (t2.counter > 1000) {
			t2.canonicalize();
		}
		T denominator;
		T ext1;
		T ext2;
		denominator = this.ring.multiply(t1.denominator, t2.denominator);
		ext1 = t2.denominator;
		ext2 = t1.denominator;
		T numerator = ring.add(ring.multiply(t1.numerator, ext1), ring.multiply(t2.numerator, ext2));
		return new Fraction<>(this.ring, numerator, denominator, t1.counter + t2.counter + 1);
	}

	@Override
	public Fraction<T> negative(Fraction<T> t) {
		return new Fraction<T>(this.ring, this.ring.negative(t.numerator), t.denominator, t.counter);
	}

	@Override
	public Fraction<T> multiply(Fraction<T> t1, Fraction<T> t2) {
		if (t1.numerator.equals(t2.denominator)) {
			return new Fraction<>(ring, t2.numerator, t1.denominator, Math.max(t1.counter, t2.counter));
		}
		if (t1.denominator.equals(t2.numerator)) {
			return new Fraction<>(ring, t1.numerator, t2.denominator, Math.max(t1.counter, t2.counter));
		}
		return new Fraction<>(ring, ring.multiply(t1.numerator, t2.numerator),
				ring.multiply(t1.denominator, t2.denominator), Math.max(t1.counter, t2.counter));
	}

	@Override
	public Fraction<T> inverse(Fraction<T> t) {
		return new Fraction<T>(this.ring, t.denominator, t.numerator, t.counter);
	}

	@Override
	public FactorizationResult<Polynomial<Fraction<T>>, Fraction<T>> factorization(
			UnivariatePolynomial<Fraction<T>> t) {
		if (!ring.isUniqueFactorizationDomain()) {
			throw new ArithmeticException("Not a UFD!");
		}
		UnivariatePolynomialRing<Fraction<T>> polynomialRing = getUnivariatePolynomialRing();
		Fraction<T> unit = t.leadingCoefficient();
		t = polynomialRing.normalize(t);
		T lcm = ring.one();
		for (int i = 0; i <= t.degree(); i++) {
			lcm = ring.lcm(lcm, t.univariateCoefficient(i).getDenominator());
		}
		t = polynomialRing.toUnivariate(polynomialRing.scalarMultiply(getEmbedding(lcm), t));
		UnivariatePolynomial<T> overRing = ring.getUnivariatePolynomialRing().getEmbedding(t, new MathMap<>() {
			@Override
			public T evaluate(Fraction<T> t) {
				return t.asInteger();
			}
		});
		FactorizationResult<Polynomial<T>, T> overRingFactors = ring.factorization(overRing);
		SortedMap<Polynomial<Fraction<T>>, Integer> factors = new TreeMap<>();
		for (Polynomial<T> overRingFactor : overRingFactors.primeFactors()) {
			int multiplicity = overRingFactors.multiplicity(overRingFactor);
			Polynomial<Fraction<T>> factor = polynomialRing.getEmbedding(overRingFactor, new MathMap<>() {
				@Override
				public Fraction<T> evaluate(T t) {
					return getEmbedding(t);
				}
			});
			factors.put(polynomialRing.normalize(factor), multiplicity);
		}
		return new FactorizationResult<>(unit, factors);
	}

	@Override
	public Fraction<T> getNumerator(Fraction<T> t) {
		return getEmbedding(t.getNumerator());
	}

	@Override
	public Fraction<T> getDenominator(Fraction<T> t) {
		return getEmbedding(t.getDenominator());
	}

	@Override
	public Fraction<T> getRandomElement() {
		T nonzero = this.ring.getRandomElement();
		while (nonzero.equals(this.ring.zero())) {
			nonzero = this.ring.getRandomElement();
		}
		return new Fraction<T>(this.ring, this.ring.getRandomElement(), nonzero, 1);
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return this.ring.getNumberOfElements();
	}

	@Override
	public Iterator<Fraction<T>> iterator() {
		throw new InfinityException();
	}

	public Fraction<T> getEmbedding(T t) {
		return new Fraction<T>(this.ring, t, this.ring.one(), 0);
	}

	public Fraction<T> getFraction(T numerator, T denominator) {
		return new Fraction<>(ring, numerator, denominator, 1);
	}

	public static class Fraction<T extends Element<T>> extends AbstractElement<Fraction<T>> {
		private Ring<T> ring;
		private T numerator;
		private T denominator;
		private boolean canonical;
		private int counter;

		private Fraction(Ring<T> ring, T numerator, T denominator, int counter) {
			if (denominator.equals(ring.zero())) {
				throw new ArithmeticException("Division by 0!");
			}
			this.ring = ring;
			if (numerator.equals(ring.zero())) {
				denominator = ring.one();
			}
			if (numerator.equals(denominator)) {
				numerator = ring.one();
				denominator = ring.one();
			}
			this.numerator = numerator;
			this.denominator = denominator;
			this.canonical = false;
			this.counter = counter;
		}

		public void canonicalize() {
			if (canonical) {
				return;
			}
			this.canonical = true;
			this.counter = 0;
			if (this.denominator.equals(ring.one())) {
				return;
			}
			if (this.ring.isUniqueFactorizationDomain()) {
				T gcd = this.ring.gcd(this.numerator, this.denominator);
				this.numerator = this.ring.divideChecked(this.numerator, gcd);
				this.denominator = this.ring.divideChecked(this.denominator, gcd);
			}
			T unit = this.ring.projectToUnit(this.denominator);
			if (!unit.equals(this.ring.one())) {
				this.numerator = this.ring.divide(this.numerator, unit);
				this.denominator = this.ring.divide(this.denominator, unit);
			}

		}

		public T getNumerator() {
			canonicalize();
			return this.numerator;
		}

		public T getDenominator() {
			canonicalize();
			return this.denominator;
		}

		public T getNumeratorDirect() {
			return this.numerator;
		}

		public T getDenominatorDirect() {
			return this.denominator;
		}

		public T asInteger() {
			QuotientAndRemainderResult<T> qr = ring.quotientAndRemainder(numerator, denominator);
			if (!qr.getRemainder().equals(ring.zero())) {
				throw new ArithmeticException("not an integer!");
			}
			return qr.getQuotient();
		}

		public String toString() {
			canonicalize();
			if (this.denominator.equals(this.ring.one())) {
				return this.numerator.toString();
			}
			String numeratorString = numerator.toString();
			boolean spacesInNumerator = numeratorString.contains(" ");
			String denominatorString = denominator.toString();
			boolean spacesInDenominator = denominatorString.contains(" ");
			return (spacesInNumerator ? "(" : "") + numeratorString + (spacesInNumerator ? ")" : "") + "/"
					+ (spacesInDenominator ? "(" : "") + denominatorString + (spacesInDenominator ? ")" : "");
		}

		@Override
		public boolean equals(Object o) {
			if (!(o instanceof Fraction<?>)) {
				return false;
			}
			@SuppressWarnings("unchecked")
			Fraction<T> other = (Fraction<T>) o;
			return ring.multiply(numerator, other.denominator).equals(ring.multiply(other.numerator, denominator));
		}

		@Override
		public int compareTo(Fraction<T> e) {
			canonicalize();
			e.canonicalize();
			return this.ring.multiply(this.numerator, e.denominator)
					.compareTo(this.ring.multiply(e.numerator, this.denominator));
		}

	}
}