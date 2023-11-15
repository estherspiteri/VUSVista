import {
  ILoadAllVusResponse,
  IStoreAndVerifyVusFileRequest,
  IStoreAndVerifyVusFileResponse,
} from "./vus.dto";

export class VusService {
  async storeAndVerifyVusFile(
    input: IStoreAndVerifyVusFileRequest
  ): Promise<IStoreAndVerifyVusFileResponse> {
    let data = new FormData();
    data.append("file", input.vusFile);

    const result: IStoreAndVerifyVusFileResponse = await fetch(`/vus/file`, {
      method: "POST",
      body: data,
      // headers: {
      // "Content-Type": "multipart/form-data",
      // },
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async loadAllVus(): Promise<ILoadAllVusResponse> {
    const result: ILoadAllVusResponse = await fetch(`/vus/view`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const vusService = new VusService();
