package fields.local;

import fields.helper.AbstractElement;
import fields.interfaces.Element;

public class Value extends AbstractElement<Value> implements Element<Value> {
	final public static Value INFINITY = new Value();
	final public static Value ZERO = new Value(0);
	final public static Value ONE = new Value(1);

	private boolean infinity;
	private int value;

	private Value() {
		infinity = true;
	}

	public Value(int value) {
		infinity = false;
		this.value = value;
	}

	public int value() {
		if (infinity) {
			throw new ArithmeticException("Infinite Value!");
		}
		return value;
	}

	public boolean isInfinite() {
		return infinity;
	}

	@Override
	public String toString() {
		if (infinity) {
			return "inf";
		}
		return Integer.toString(value);
	}

	@Override
	public int compareTo(Value o) {
		if (this.infinity && o.infinity) {
			return 0;
		}
		if (this.infinity) {
			return 1;
		}
		if (o.infinity) {
			return -1;
		}
		return this.value - o.value;
	}
	
	public Value min(Value o) {
		if (this.compareTo(o) < 0) {
			return this;
		}
		return o;
	}
	
	public Value max(Value o) {
		if (this.compareTo(o) < 0) {
			return o;
		}
		return this;
	}
}
