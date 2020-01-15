package fields;

import fields.FieldOfFractions.Fraction;

public class FieldOfFractions<T extends Element> extends AbstractField<Fraction<T>> implements Field<Fraction<T>> {
	private Ring<T> ring;
	public FieldOfFractions(Ring<T> ring) {
		this.ring = ring;
		if (!this.ring.isIntegral())
			throw new RuntimeException("Object does not exist");
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
	public int characteristic() {
		return this.ring.characteristic();
	}
	@Override
	public Fraction<T> add(Fraction<T> t1, Fraction<T> t2) {
		T numerator = this.ring.add(this.ring.multiply(t1.getNumerator(), t2.getDenominator()), this.ring.multiply(t2.getNumerator(), t1.getDenominator()));
		T denominator = this.ring.multiply(t1.getDenominator(), t2.getDenominator());
		return new Fraction<T>(this.ring, numerator, denominator);
	}
	@Override
	public Fraction<T> negative(Fraction<T> t) {
		return new Fraction<T>(this.ring, this.ring.negative(t.getNumerator()), t.getDenominator());
	}
	@Override
	public Fraction<T> multiply(Fraction<T> t1, Fraction<T> t2) {
		T numerator = this.ring.multiply(t1.getNumerator(), t2.getNumerator());
		T denominator = this.ring.multiply(t1.getDenominator(), t2.getDenominator());
		return new Fraction<T>(this.ring, numerator, denominator);
	}
	@Override
	public Fraction<T> inverse(Fraction<T> t) {
		return new Fraction<T>(this.ring, t.getDenominator(), t.getNumerator());
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
	public int getNumberOfElements() throws InfinityException {
		return this.ring.getNumberOfElements();
	}
	@Override
	public Iterable<Fraction<T>> getElements() throws InfinityException {
		throw new InfinityException();
	}
	public Fraction<T> getEmbedding(T t) {
		return new Fraction<T>(this.ring, t, this.ring.one());
	}
	public static class Fraction<T extends Element> implements Element {
		private Ring<T> ring;
		private T numerator;
		private T denominator;
		
		private Fraction(Ring<T> ring, T numerator, T denominator) {
			this.ring = ring;
			this.numerator = numerator;
			this.denominator = denominator;
			if (this.denominator.compareTo(this.ring.zero()) < 0) {
				this.numerator = this.ring.negative(this.numerator);
				this.denominator = this.ring.negative(this.denominator);
			}
			if (this.ring.isEuclidean()) {
				T gcd = this.ring.gcd(this.numerator, this.denominator);
				this.numerator = this.ring.quotientAndRemainder(this.numerator, gcd).get(0);
				this.denominator = this.ring.quotientAndRemainder(this.denominator, gcd).get(0);
			}
		}
		
		public T getNumerator() {
			return this.numerator;
		}

		public T getDenominator() {
			return this.denominator;
		}

		public boolean equals(Object o) {
			if (!(o instanceof FieldOfFractions.Fraction))
				return false;
			@SuppressWarnings("unchecked")
			Fraction<T> e = (Fraction<T>)o;
			return this.ring.multiply(this.numerator,  e.denominator).equals(
					this.ring.multiply(e.numerator,  this.denominator));
		}
		public String toString() {
			if (this.denominator.equals(this.ring.one()))
				return this.numerator.toString();
			return this.numerator.toString() + "/" + this.denominator.toString();
		}

		@Override
		public int compareTo(Element o) {
			if (!(o instanceof FieldOfFractions.Fraction))
				throw new RuntimeException();
			@SuppressWarnings("unchecked")
			Fraction<T> e = (Fraction<T>)o;
			return this.ring.multiply(this.numerator,  e.denominator).compareTo(
					this.ring.multiply(e.numerator,  this.denominator));
		}
	}
}