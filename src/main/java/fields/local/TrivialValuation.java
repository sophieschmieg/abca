package fields.local;

import fields.interfaces.Element;
import fields.interfaces.DiscreteValuationField.Valuation;
import fields.interfaces.Ring;

public class TrivialValuation<T extends Element<T>> implements Valuation<T> {
	private Ring<T> ring;

	public TrivialValuation(Ring<T> ring) {
		this.ring = ring;
	}

	@Override
	public Value valuation(T t) {
		return t.equals(ring.zero()) ? Value.INFINITY : Value.ZERO;
	}

}
