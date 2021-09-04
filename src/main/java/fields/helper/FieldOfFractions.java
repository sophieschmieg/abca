package fields.helper;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.helper.FieldOfFractions.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ring;

public class FieldOfFractions<T extends Element<T>> extends AbstractField<Fraction<T>> implements Field<Fraction<T>> {
	private Ring<T> ring;

	public FieldOfFractions(Ring<T> ring) {
		this.ring = ring;
		if (!this.ring.isIntegral())
			throw new RuntimeException("Object does not exist");
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
			return new Fraction<>(this.ring, this.ring.add(t1.numerator, t2.numerator), t1.denominator);
		}
		T denominator;
		T ext1;
		T ext2;
//		if (ring.isUniqueFactorizationDomain()) {
//			denominator = ring.lcm(t1.getDenominator(), t2.getDenominator());
//			ext1 = ring.divideChecked(denominator, t1.getDenominator());
//			ext2 = ring.divideChecked(denominator, t2.getDenominator());
//		} else {
			denominator = this.ring.multiply(t1.denominator, t2.denominator);
			ext1 = t2.denominator;
			ext2 = t1.denominator;
		//}
		T numerator = ring.add(ring.multiply(t1.numerator, ext1), ring.multiply(t2.numerator, ext2));
		return new Fraction<>(this.ring, numerator, denominator);
	}

	@Override
	public Fraction<T> negative(Fraction<T> t) {
		return new Fraction<T>(this.ring, this.ring.negative(t.numerator), t.denominator);
	}

	@Override
	public Fraction<T> multiply(Fraction<T> t1, Fraction<T> t2) {
		if (t1.numerator.equals(t2.denominator)) {
			return new Fraction<>(ring, t2.numerator, t1.denominator);
		}
		if (t1.denominator.equals(t2.numerator)) {
			return new Fraction<>(ring, t1.numerator, t2.denominator);
		}
		return new Fraction<>(ring, ring.multiply(t1.numerator, t2.numerator),
				ring.multiply(t1.denominator, t2.denominator));
	}

	@Override
	public Fraction<T> inverse(Fraction<T> t) {
		return new Fraction<T>(this.ring, t.denominator, t.numerator);
	}

	@Override
	public Fraction<T> numerator(Fraction<T> t) {
		return getEmbedding(t.getNumerator());
	}

	@Override
	public Fraction<T> denominator(Fraction<T> t) {
		return getEmbedding(t.getDenominator());
	}

	@Override
	public Fraction<T> getRandomElement() {
		T nonzero = this.ring.getRandomElement();
		while (nonzero.equals(this.ring.zero())) {
			nonzero = this.ring.getRandomElement();
		}
		return new Fraction<T>(this.ring, this.ring.getRandomElement(), nonzero);
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
		return new Fraction<T>(this.ring, t, this.ring.one());
	}

	public Fraction<T> getFraction(T numerator, T denominator) {
		return new Fraction<>(ring, numerator, denominator);
	}

	public static class Fraction<T extends Element<T>> extends AbstractElement<Fraction<T>> {
		private Ring<T> ring;
		private T numerator;
		private T denominator;
		private boolean canonical;

		private Fraction(Ring<T> ring, T numerator, T denominator) {
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
		}

		public void canonicalize() {
			if (canonical) {
				return;
			}
			this.canonical = true;
			if (this.denominator.equals(ring.one())) {
				return;
			}
			if (this.ring.isUniqueFactorizationDomain()) {
				T gcd = this.ring.gcd(this.numerator, this.denominator);
				this.numerator = this.ring.divideChecked(this.numerator, gcd);
				this.denominator = this.ring.divideChecked(this.denominator, gcd);
			}
			if (this.denominator.compareTo(this.ring.zero()) < 0) {
				this.numerator = this.ring.negative(this.numerator);
				this.denominator = this.ring.negative(this.denominator);
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