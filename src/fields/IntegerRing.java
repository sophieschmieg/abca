package fields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import fields.IntegerRing.IntegerRingElement;

public class IntegerRing extends AbstractRing<IntegerRingElement> {

	public static class IntegerRingElement implements Element {
		private BigInteger value;

		public IntegerRingElement(int value) {
			this.value = BigInteger.valueOf(value);
		}
		
		public IntegerRingElement(BigInteger value) {
			this.value = value;
		}

		@Override
		public int compareTo(Element other) {
			return this.value.compareTo(((IntegerRingElement) other).value);
		}

		public boolean equals(Object o) {
			return this.value.equals(((IntegerRingElement) o).value);
		}

		public BigInteger getValue() {
			return value;
		}

		public String toString() {
			return this.value.toString();
		}
	}

	@Override
	public IntegerRingElement getRandomElement() {
		return new IntegerRingElement((int)Math.round(new Random().nextGaussian() * 10.0));
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
	public Iterator<IntegerRingElement> iterator() {
		throw new InfinityException();
	}

	@Override
	public IntegerRingElement zero() {
		return new IntegerRingElement(0);
	}

	@Override
	public IntegerRingElement one() {
		return new IntegerRingElement(1);
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public IntegerRingElement add(IntegerRingElement t1, IntegerRingElement t2) {
		return new IntegerRingElement(t1.value.add(t2.value));
	}

	@Override
	public IntegerRingElement negative(IntegerRingElement t) {
		return new IntegerRingElement(t.value.negate());
	}

	@Override
	public IntegerRingElement multiply(IntegerRingElement t1, IntegerRingElement t2) {
		return new IntegerRingElement(t1.value.multiply(t2.value));
	}
	
	@Override
	public IntegerRingElement getInteger(BigInteger t) {
		return new IntegerRingElement(t);
	}

	@Override
	public boolean isUnit(IntegerRingElement t) {
		return t.value.equals(BigInteger.ONE) || t.value.equals(BigInteger.ONE.negate());
	}

	@Override
	public BigInteger getNumberOfUnits() {
		return BigInteger.TWO;
	}

	@Override
	public Iterable<IntegerRingElement> getUnits() {
		List<IntegerRingElement> list = new ArrayList<>();
		list.add(new IntegerRingElement(1));
		list.add(new IntegerRingElement(-1));
		return list;
	}

	@Override
	public IntegerRingElement inverse(IntegerRingElement t) {
		return t;
	}

	@Override
	public boolean isIntegral() {
		return true;
	}
	
	@Override
	public boolean isZeroDivisor(IntegerRingElement t) {
		return t.value.equals(BigInteger.ZERO);
	}

	@Override
	public boolean isEuclidean() {
		return true;
	}
	
	@Override
	public boolean isDivisible(IntegerRingElement dividend,
			IntegerRingElement divisor) {
		return dividend.value.mod(divisor.value).equals(BigInteger.ZERO);
	}
	
	@Override
	public BigInteger euclidMeasure(IntegerRingElement t) {
		return t.value.abs();
	}

	@Override
	public List<IntegerRingElement> quotientAndRemainder(IntegerRingElement dividend,
			IntegerRingElement divisor) {
		BigInteger[] result = dividend.value.divideAndRemainder(divisor.value);
		List<IntegerRingElement> list = new ArrayList<>();
		list.add(new IntegerRingElement(result[0]));
		list.add(new IntegerRingElement(result[1]));
		return list;
	}
}
