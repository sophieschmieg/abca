package util;

import java.util.ArrayList;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.MathMap;
import fields.vectors.Vector;

public class ProductVectorMap<T extends Element<T>, S extends Element<S>> implements MathMap<Vector<T>, Vector<S>> {
private MathMap<T, S> map;
public ProductVectorMap(MathMap<T, S> map) {
	this.map = map;
}
	
	@Override
	public Vector<S> evaluate(Vector<T> t) {
		List<S> results = new ArrayList<>();
		for (T e : t.asList()) {
			results.add(map.evaluate(e));
		}
		return new Vector<>(results);
	}

}
