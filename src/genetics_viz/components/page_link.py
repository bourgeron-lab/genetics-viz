"""Shareable page state: the open tab and the filters, kept in the URL's query.

A curator sends the address bar to a colleague, who opens the same tab with
the same filters. Pages read their state from the query string they were
opened with, and write it back with ``history.replaceState`` as the user
changes it, so the address bar is always the link to the current view and
the browser history gets no entry per click.

Only values that differ from a page's defaults are written, so a page opened
with no query is in its usual state, and a link stays short. Lists are
comma-joined: none of the values linked (sample IDs, geneset keys, VEP terms,
statuses, call types) contains a comma.

The query reaches the server, and so its access log, as the path already
does. That is the price of reading the state as the page is built.
"""

import json
from typing import Any, Dict, Iterable, List, Mapping, Optional
from urllib.parse import urlencode

from genetics_viz.utils.vep import VEP_CONSEQUENCES


class PageLink:
    """A page's query-string state, read at build time and kept in the URL."""

    def __init__(self, path: str, params: Mapping[str, str]) -> None:
        self.path = path
        self.params: Dict[str, str] = {k: v for k, v in params.items() if v != ""}
        self._client: Any = None

    @classmethod
    def from_request(cls, request: Any) -> "PageLink":
        """The state of the page being built (``request`` from the page function)."""
        from nicegui import context

        link = cls(request.url.path, dict(request.query_params))
        link._client = context.client
        return link

    # ---- reading ----

    def get(self, key: str, default: Optional[str] = None) -> Optional[str]:
        return self.params.get(key, default)

    def get_list(self, key: str) -> Optional[List[str]]:
        """The comma-joined list, or ``None`` when the link does not set it.

        ``-`` stands for an empty list, which a link must tell apart from
        "not set" (an emptied filter differs from the default one).
        """
        raw = self.params.get(key)
        if raw is None:
            return None
        if raw == "-":
            return []
        return [v for v in raw.split(",") if v]

    def get_bool(self, key: str, default: bool) -> bool:
        raw = self.params.get(key)
        if raw is None:
            return default
        return raw not in ("0", "false", "no")

    # ---- writing ----

    def set(self, key: str, value: Any, default: Any = None) -> None:
        """Record ``value``; drop the key when it equals ``default`` (or is None).

        Lists compare as sets with their default: a selection's order is not
        part of the view.
        """
        if value is None or _same(value, default):
            self.params.pop(key, None)
        elif isinstance(value, bool):
            self.params[key] = "1" if value else "0"
        elif isinstance(value, (list, tuple, set)):
            self.params[key] = ",".join(str(v) for v in value) or "-"
        else:
            self.params[key] = str(value)

    def drop(self, *keys: str) -> None:
        for key in keys:
            self.params.pop(key, None)

    def url(self) -> str:
        """Path and query of the current state."""
        if not self.params:
            return self.path
        return f"{self.path}?{urlencode(self.params, safe=',:')}"

    def push(self) -> None:
        """Put :meth:`url` in the address bar, replacing the history entry."""
        if self._client is None:
            return
        from nicegui import ui

        with self._client:
            ui.run_javascript(
                f"history.replaceState(history.state, '', {_js(self.url())})"
            )


def _same(value: Any, default: Any) -> bool:
    if isinstance(value, (list, tuple, set)) and isinstance(
        default, (list, tuple, set)
    ):
        return set(value) == set(default)
    return value == default


def _js(text: str) -> str:
    return json.dumps(text)


#: VEP impact levels, the shorthand of a consequence selection.
_LEVELS = ("HIGH", "MODERATE", "LOW", "MODIFIER")


def encode_consequences(selected: Iterable[str]) -> List[str]:
    """A consequence selection, as whole impact levels where it can be.

    ``HIGH`` stands for every HIGH term; terms of a partly selected level are
    listed one by one. Keeps a link with the usual presets short.
    """
    chosen = set(selected)
    out: List[str] = []
    for level in _LEVELS:
        terms = {t for t, (imp, _) in VEP_CONSEQUENCES.items() if imp == level}
        if terms and terms <= chosen:
            out.append(level)
            chosen -= terms
    return out + sorted(t for t in chosen if t in VEP_CONSEQUENCES)


def decode_consequences(values: Iterable[str]) -> List[str]:
    """The terms of :func:`encode_consequences`'s output; unknown ones dropped."""
    out: List[str] = []
    for v in values:
        if v in _LEVELS:
            out.extend(t for t, (imp, _) in VEP_CONSEQUENCES.items() if imp == v)
        elif v in VEP_CONSEQUENCES:
            out.append(v)
    return list(dict.fromkeys(out))


def encode_sorting(sorting: List[Dict[str, Any]]) -> Optional[str]:
    """``col`` or ``-col`` (descending), comma-joined; None when unsorted."""
    if not sorting:
        return None
    return ",".join(("-" if s.get("desc") else "") + str(s["id"]) for s in sorting)


def decode_sorting(raw: Optional[str]) -> List[Dict[str, Any]]:
    if not raw:
        return []
    out = []
    for part in raw.split(","):
        if part:
            desc = part.startswith("-")
            out.append({"id": part.lstrip("-"), "desc": desc})
    return out


def safe_next(target: Optional[str]) -> Optional[str]:
    """A login ``next`` target, only when it is a path on this site.

    ``//host`` and ``/\\host`` are other sites to a browser, and a full URL
    is never accepted, so the login cannot be used to send someone elsewhere.
    """
    if not target or not target.startswith("/"):
        return None
    if target.startswith("//") or target.startswith("/\\"):
        return None
    if target.split("?", 1)[0] == "/login":
        return None
    return target


def link_members(
    link: PageLink,
    selected: Dict[str, List[str]],
    member_ids: List[str],
    refreshers: List[Any],
) -> None:
    """Link a family page's member selection (``m``), before the selector is built.

    The link's members that are in the family are selected, in pedigree
    order; a link naming none of them keeps the default (everyone). The link
    follows every selection change through the page's refreshers.
    """
    linked = link.get_list("m")
    if linked is not None:
        chosen = [m for m in member_ids if m in linked]
        if chosen:
            selected["value"] = chosen

    def write() -> None:
        link.set("m", selected["value"], default=member_ids)
        link.push()

    refreshers.append(write)
