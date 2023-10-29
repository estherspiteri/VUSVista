class LitvarResponse:
  def __init__(self, data, status: int, error_msg: str | None = None,):
    self.data = data
    self.status = status
    self.error_msg = error_msg
