package fields.vectors.pivot;

import java.util.List;

import fields.interfaces.Element;

public interface PivotStrategy<T extends Element<T>> {
	public int pivot(List<T> elements);
}
