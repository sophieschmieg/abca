package fields.vectors.pivot;

import java.util.List;

import fields.floatingpoint.Reals.Real;
import fields.interfaces.Element;
import fields.interfaces.ValueField;

public class ValuePivotStrategy<T extends Element<T>> implements PivotStrategy<T> {
	private ValueField<T> field;

	public ValuePivotStrategy(ValueField<T> field) {
		this.field = field;
	}

	@Override
	public int pivot(List<T> elements) {
		Real min = null;
		int minPos = -1;
		for (int i = 0; i < elements.size(); i++) {
			if (elements.get(i).equals(field.zero())) {
				continue;
			}
			Real value = field.value(elements.get(i));
			if (min == null || value.compareTo(min) < 0) {
				minPos = i;
				min = value;
			}
		}
		return minPos;
	}
}
