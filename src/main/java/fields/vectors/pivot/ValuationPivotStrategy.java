package fields.vectors.pivot;

import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.LocalField.Valuation;
import fields.local.Value;

public class ValuationPivotStrategy<T extends Element<T>> implements PivotStrategy<T> {
	private Valuation<T> valuation;

	public ValuationPivotStrategy(Valuation<T> valuation) {
		this.valuation = valuation;
	}

	@Override
	public int pivot(List<T> elements) {
		Value minValue = Value.INFINITY;
		int pivot = -1;
		for (int i = 0; i < elements.size(); i++) {
			Value value = valuation.valuation(elements.get(i));
			if (value.compareTo(minValue) < 0) {
				pivot = i;
				minValue = value;
			}
		}
		return pivot;
	}
}
