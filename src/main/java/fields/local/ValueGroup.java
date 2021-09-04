package fields.local;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.interfaces.Monoid;
import fields.interfaces.WellOrder;

public class ValueGroup implements Monoid<Value>, WellOrder<Value> {
private static ValueGroup g = new ValueGroup();
	private  ValueGroup() {
	}
	
	public static ValueGroup g() {
		return g;
	}
	
	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public Value getRandomElement() {
		// TODO Auto-generated method stub
		return null;
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
	public Iterator<Value> iterator() {
		return new Iterator<Value>() {
			private Value i = least();

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public Value next() {
				Value next = i;
				i = new Value(i.value() + 1);
				return next;
			}
		};
	}

	@Override
	public Value least() {
		return new Value(0);
	}

	@Override
	public int compare(Value t1, Value t2) {
		return t1.compareTo(t2);
	}

	@Override
	public Value neutral() {
		return new Value(0);
	}

	@Override
	public Value operate(Value t1, Value t2) {
		if (t1.isInfinite() || t2.isInfinite()) {
			return Value.INFINITY;
		}
		return new Value(t1.value() + t2.value());
	}

	public Value multiply(int n, Value t) {
		if (t.isInfinite()) {
			return t;
		}
		return new Value(n * t.value());
	}
	
	public Value add(Value t1, Value t2) {
		return operate(t1, t2);
	}
	
	public Value subtract(Value t1, Value t2) {
		if (t1.isInfinite()) {
			return t1;
		}
		if (t2.isInfinite()) {
			throw new ArithmeticException("Cannot subtract infinity");
		}
		return new Value(t1.value() - t2.value());
	}
}
