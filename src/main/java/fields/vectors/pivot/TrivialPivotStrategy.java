package fields.vectors.pivot;

import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Ring;

public class TrivialPivotStrategy<T extends Element<T>> implements PivotStrategy<T> {
private Ring<T> ring;

public TrivialPivotStrategy(Ring<T> ring) {
	this.ring = ring;
}

@Override
public int pivot(List<T> elements) {
	for (int i = 0; i < elements.size(); i++) {
		if (!elements.get(i).equals(ring.zero())) {
			return i;
		}
	}
	return -1;
}
}
