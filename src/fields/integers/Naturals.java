package fields.integers;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.Random;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.integers.Naturals.Natural;
import fields.interfaces.WellOrder;

public class Naturals implements WellOrder<Natural> {
	private static Naturals n = new Naturals();

	public static class Natural extends AbstractElement<Natural> {
		private BigInteger value;

		public Natural(int value) {
			this(BigInteger.valueOf(value));
		}

		public Natural(BigInteger value) {
			if (value.compareTo(BigInteger.ZERO) < 0) {
				throw new ArithmeticException("Not a natural number");
			}
			this.value = value;
		}

		@Override
		public int compareTo(Natural o) {
			return value.compareTo(o.value);
		}

		@Override
		public String toString() {
			return value.toString();
		}
	}

	private Naturals() {
	}

	public static Naturals n() {
		return n;
	}

	@Override
	public Natural getRandomElement() {
		return new Natural((int) Math.abs(Math.round(new Random().nextGaussian() * 10.0)));
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
	public Iterator<Natural> iterator() {
		return new Iterator<Natural>() {
			BigInteger i = BigInteger.ZERO;

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public Natural next() {
				BigInteger next = i;
				i = i.add(BigInteger.ONE);
				return new Natural(next);
			}
			
		};
	}

	@Override
	public Natural least() {
		return new Natural(0);
	}
	
	@Override
	public Natural increment(Natural t) {
		return new Natural(t.value.add(BigInteger.ONE));
	}
}
