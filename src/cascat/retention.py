"""Independent per-stage retention, with percentages rounded up."""

from decimal import Decimal, InvalidOperation, ROUND_CEILING


def parse_retention(value):
    if value.lower() == "all":
        return None
    if value.endswith("%"):
        try:
            percent = Decimal(value[:-1])
        except InvalidOperation:
            raise ValueError("percentage must be a number in (0, 100]") from None
        if not percent.is_finite() or not 0 < percent <= 100:
            raise ValueError("percentage must be in (0, 100]")
        return f"{percent.normalize():f}%"
    number = int(value)
    if number < 1:
        raise ValueError("retention must be a positive integer, percentage, or all")
    return number


def retention_count(available, requested):
    if requested is None:
        return available
    if isinstance(requested, str):
        fraction = Decimal(requested[:-1]) / 100
        return int((available * fraction).to_integral_value(rounding=ROUND_CEILING))
    return min(available, requested)


def label(requested):
    return "all" if requested is None else requested
