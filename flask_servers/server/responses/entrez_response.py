class EntrezResponse:
  def __init__(self, data, status: int, error_msg: str):
    self.data = data
    self.status = status
    self.error_msg = error_msg
