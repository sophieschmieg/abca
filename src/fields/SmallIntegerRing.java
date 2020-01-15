package fields;

import java.util.ArrayList;
import java.util.List;

import fields.SmallIntegerRing.SmallIntegerRingElement;

public class SmallIntegerRing extends AbstractRing<SmallIntegerRingElement> {

	public static class SmallIntegerRingElement implements Element {
		private int value;

		public SmallIntegerRingElement(int value) {
			this.value = value;
		}

		@Override
		public int compareTo(Element other) {
			return this.value - ((SmallIntegerRingElement)other).value;
		}
		
		public boolean equals(Object o) {
			return this.value == ((SmallIntegerRingElement)o).value;
		}

		public int getValue() {
			return value;
		}
		
		public String toString() {
			return Integer.toString(this.value);
		}
	}

	@Override
	public SmallIntegerRingElement getRandomElement() {
		return new SmallIntegerRingElement(0);
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public int getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterable<SmallIntegerRingElement> getElements()
			throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public SmallIntegerRingElement zero() {
		return new SmallIntegerRingElement(0);
	}

	@Override
	public SmallIntegerRingElement one() {
		return new SmallIntegerRingElement(1);
	}

	@Override
	public int characteristic() {
		return 0;
	}

	@Override
	public SmallIntegerRingElement add(
			SmallIntegerRingElement t1, SmallIntegerRingElement t2) {
		return new SmallIntegerRingElement(t1.value + t2.value);
	}

	@Override
	public SmallIntegerRingElement negative(
			SmallIntegerRingElement t) {
		return new SmallIntegerRingElement(-t.value);
	}

	@Override
	public SmallIntegerRingElement multiply(
			SmallIntegerRingElement t1, SmallIntegerRingElement t2) {
		return new SmallIntegerRingElement(t1.value * t2.value);
	}

	@Override
	public boolean isUnit(SmallIntegerRingElement t) {
		return t.value == 1 || t.value == -1;
	}

	@Override
	public int getNumberOfUnits() {
		return 2;
	}

	@Override
	public Iterable<SmallIntegerRingElement> getUnits() {
		List<SmallIntegerRingElement> list = new ArrayList<SmallIntegerRing.SmallIntegerRingElement>();
		list.add(new SmallIntegerRingElement(1));
		list.add(new SmallIntegerRingElement(-1));
		return list;
	}

	@Override
	public SmallIntegerRingElement inverse(
			SmallIntegerRingElement t) {
		return t;
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isEuclidean() {
		return true;
	}

	@Override
	public List<SmallIntegerRingElement> quotientAndRemainder(
			SmallIntegerRingElement dividend,
			SmallIntegerRingElement divisor) {
		int q = dividend.value / divisor.value;
		int r = dividend.value - q * divisor.value;
		List<SmallIntegerRingElement> list = new ArrayList<SmallIntegerRing.SmallIntegerRingElement>();
		list.add(new SmallIntegerRingElement(q));
		list.add(new SmallIntegerRingElement(r));
		return list;
	}
}
