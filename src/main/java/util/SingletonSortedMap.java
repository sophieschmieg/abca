package util;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;

public class SingletonSortedMap<K, V> implements SortedMap<K, V> {
	private K key;
	private V value;
	private Comparator<? super K> comparator;

	public static <K extends Comparable<K>, V> SingletonSortedMap<K, V> map(K key, V value) {
		return new SingletonSortedMap<>(key, value, Comparator.naturalOrder());
	}

	public static <K, V> SingletonSortedMap<K, V> map(K key, V value, Comparator<? super K> comparator) {
		return new SingletonSortedMap<>(key, value, comparator);
	}

	private SingletonSortedMap(K key, V value, Comparator<? super K> comparator) {
		this.key = key;
		this.value = value;
		this.comparator = comparator;
	}

	@Override
	public int size() {
		return 1;
	}

	@Override
	public boolean isEmpty() {
		return false;
	}

	@Override
	public boolean containsKey(Object key) {
		return key.equals(this.key);
	}

	@Override
	public boolean containsValue(Object value) {
		return value.equals(value);
	}

	@Override
	public V get(Object key) {
		if (key.equals(this.key)) {
			return value;
		}
		return null;
	}

	@Override
	public V put(K key, V value) {
if (key.equals(this.key)) {
	this.value = value;
}
throw new UnsupportedOperationException("Singleton Map cannot add items, only change");
	}

	@Override
	public V remove(Object key) {
		throw new UnsupportedOperationException("Singleton Map cannot remove items, only change");
	}

	@Override
	public void putAll(Map<? extends K, ? extends V> m) {
	for (K key : m.keySet()) {
		put(key, m.get(key));
	}
	}

	@Override
	public void clear() {
		throw new UnsupportedOperationException("Singleton Map cannot remove items, only change");
		
	}

	@Override
	public Comparator<? super K> comparator() {
		return comparator;
	}

	@Override
	public SortedMap<K, V> subMap(K fromKey, K toKey) {
		if (comparator.compare(key, fromKey) >= 0 && comparator.compare(key, toKey) < 0) {
			return this;
		}
		return Collections.emptySortedMap();
	}

	@Override
	public SortedMap<K, V> headMap(K toKey) {
		if ( comparator.compare(key, toKey) < 0) {
			return this;
		}
		return Collections.emptySortedMap();
	}

	@Override
	public SortedMap<K, V> tailMap(K fromKey) {
		if (comparator.compare(key, fromKey) >= 0) {
			return this;
		}
		return Collections.emptySortedMap();
	}

	@Override
	public K firstKey() {
		return key;
	}

	@Override
	public K lastKey() {
		return key;
	}

	@Override
	public Set<K> keySet() {
		return Collections.singleton(key);
	}

	@Override
	public Collection<V> values() {
		return Collections.singleton(value);
		}

	@Override
	public Set<Entry<K, V>> entrySet() {
		return Collections.singleton(new Entry<>() {

			@Override
			public K getKey() {
				return key;
			}

			@Override
			public V getValue() {
				return value;
			}

			@Override
			public V setValue(V value) {
				return SingletonSortedMap.this.value = value;
			}});
	}
}
