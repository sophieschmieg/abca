package fields.vectors.pivot;

import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Ring;

public class DivisionPivotStrategy<T extends Element<T>> implements PivotStrategy<T> {
	private Ring<T> ring;

	public DivisionPivotStrategy(Ring<T> ring) {
		this.ring = ring;
	}

	@Override
	public int pivot(List<T> elements) {
		T minElement = null;
		int firstNonZero = -1;
		int pivot = -1;
		for (int i = 0; i < elements.size(); i++) {
			T element = elements.get(i);
			if (element.equals(ring.zero())) {
				continue;
			}
			if (minElement == null) {
				firstNonZero = i;
				pivot = i;
				minElement = element;
			} else if (ring.isDivisible(minElement, element)) {
				pivot = i;
				minElement = element;
			} else if (!ring.isDivisible(element, minElement)) {
				return firstNonZero;
			}
		}
		return pivot;
	}

}
