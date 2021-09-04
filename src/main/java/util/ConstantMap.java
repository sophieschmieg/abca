package util;

import fields.interfaces.Element;
import fields.interfaces.MathMap;

public class ConstantMap<T extends Element<T>, S extends Element<S>> implements MathMap<T, S> {
	private S result;

	public ConstantMap(S result) {
		this.result = result;
	}

	@Override
	public S evaluate(T t) {
		return result;
	}

}
