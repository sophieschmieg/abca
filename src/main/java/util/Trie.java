package util;

import java.util.Map;
import java.util.Optional;
import java.util.TreeMap;

public class Trie {
	public class Node {
		private Character ch;
		private boolean valid;
		private Node parent;
		private Map<Character, Node> children;

		private Node() {
			this.parent = null;
			this.valid = false;
			this.children = new TreeMap<>();
		}

		private Node(Node parent, Character ch) {
			this.ch = ch;
			this.valid = false;
			this.parent = parent;
			this.children = new TreeMap<>();
		}

		private void insert(String in) {
			if (in.length() == 0) {
				valid = true;
				if (children.size() > 0) {
					prefixFree = false;
				}
				return;
			}
			char head = in.charAt(0);
			String tail = in.substring(1);
			Node child = children.get(head);
			if (child == null) {
				child = new Node(this, head);
				children.put(head, child);
			}
			child.insert(tail);
		}

		public Node child(char character) {
			return children.get(character);
		}

		public boolean isValid() {
			return valid;
		}

		public String reconstruct() {
			if (!valid) {
				throw new RuntimeException("Invalid node!");
			}
			StringBuilder build = new StringBuilder();
			Node node = this;
			while (node.parent != null) {
				build.append(node.ch);
				node = node.parent;
			}
			build.reverse();
			return build.toString();
		}
	}

	private Node root;
	private boolean prefixFree;

	public Trie() {
		this.root = new Node();
		this.prefixFree = true;
	}

	public Trie(String[] elements) {
		this();
		for (String element : elements) {
			insert(element);
		}
	}

	public void insert(String text) {
		root.insert(text);
	}
	
	public Node getRootNode() {
		return root;
	}

	public Optional<Node> fetch(String text) {
		Node node = root;
		for (int i = 0; i < text.length(); i++) {
			node = node.child(text.charAt(i));
			if (node == null) {
				return Optional.empty();
			}
		}
		return Optional.of(node);
	}

	public boolean contains(String text) {
		Optional<Node> result = fetch(text);
		return result.isPresent() && result.get().isValid();
	}

	public boolean containsPrefix(String prefix) {
		return fetch(prefix).isPresent();
	}

	public boolean isPrefixFree() {
		return prefixFree;
	}
}
