package util;

import fields.interfaces.Element;
import fields.interfaces.MathMap;

public class Identity<T extends Element<T>> implements MathMap<T, T> {

	@Override
	public T evaluate(T t) {
		return t;
	}

}
