import {
  IGetVusRequest,
  IGetVusResponse,
  ILoadAllVusResponse,
  IStoreAndVerifyVusFileRequest,
  IStoreAndVerifyVusFileResponse,
  IUploadVusRequest,
  IUploadVusResponse,
  IVerifyGeneRequest,
  IVerifyGeneResponse,
} from "./vus.dto";

export class VusService {
  async uploadVus(
    input: IUploadVusRequest
  ): Promise<IUploadVusResponse> {
    let data = new FormData();

    const vusJsonData = input.vus
      ? JSON.stringify(input.vus)
      : "";

    // Append the JSON string as a blob to the FormData
    data.append("vus", vusJsonData);

    const result: IUploadVusResponse = await fetch(`/vus/upload`, {
      method: "POST",
      body: data,
      // headers: {
      // "Content-Type": "multipart/form-data",
      // },
      cache: "no-store", //TODO: is it needed?
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async storeAndVerifyVusFile(
    input: IStoreAndVerifyVusFileRequest
  ): Promise<IStoreAndVerifyVusFileResponse> {
    let data = new FormData();
    data.append("file", input.vusFile);

    const multipleGenesSelectionJsonData = input.multipleGenesSelection
      ? JSON.stringify(input.multipleGenesSelection)
      : "";

    // Append the JSON string as a blob to the FormData
    data.append("multipleGenesSelection", multipleGenesSelectionJsonData);

    const result: IStoreAndVerifyVusFileResponse = await fetch(`/vus/file`, {
      method: "POST",
      body: data,
      // headers: {
      // "Content-Type": "multipart/form-data",
      // },
      cache: "no-store", //TODO: is it needed?
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

  async getVus(input: IGetVusRequest): Promise<IGetVusResponse> {
    const result: IGetVusResponse = await fetch(`/vus/view/${input.vusId}`, {
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

  async verifyGene(input: IVerifyGeneRequest): Promise<IVerifyGeneResponse> {
    const result: IVerifyGeneResponse = await fetch(
      `/vus/gene/${input.geneName}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const vusService = new VusService();
