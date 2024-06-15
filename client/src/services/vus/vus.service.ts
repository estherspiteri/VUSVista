import {
  IAddNewSampleRequest,
  IAddNewSampleResponse,
  IAddPublicationsRequest,
  IAddPublicationsResponse,
  IAddSamplesRequest,
  IAddSamplesResponse,
  IDeleteVariantRequest,
  IDeleteVariantResponse,
  IGetAllAcmgRulesResponse,
  IGetClinvarUpdatesRequest,
  IGetClinvarUpdatesResponse,
  IGetPublicationUpdatesRequest,
  IGetPublicationUpdatesResponse,
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
  async uploadVus(input: IUploadVusRequest): Promise<IUploadVusResponse> {
    let data = new FormData();

    const vusJsonData = input.vus ? JSON.stringify(input.vus) : "";

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
      cache: "no-cache",
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

  async getAllAcmgRules(): Promise<IGetAllAcmgRulesResponse> {
    const result: IGetAllAcmgRulesResponse = await fetch(
      `/vus/all-acmg-rules`,
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

  async getClinvarUpdates(
    input: IGetClinvarUpdatesRequest
  ): Promise<IGetClinvarUpdatesResponse> {
    const result: IGetClinvarUpdatesResponse = await fetch(
      `/vus/get_clinvar_updates/${input.clinvarId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
        cache: "no-cache",
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getPublicationUpdates(
    input: IGetPublicationUpdatesRequest
  ): Promise<IGetPublicationUpdatesResponse> {
    const result: IGetPublicationUpdatesResponse = await fetch(
      `/vus/get_publication_updates/${input.variantId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
        cache: "no-cache",
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async addPublications(
    input: IAddPublicationsRequest
  ): Promise<IAddPublicationsResponse> {
    let data = new FormData();

    const publicationsJsonData = input.publicationUrls
      ? JSON.stringify(input.publicationUrls)
      : "";

    // Append the JSON string as a blob to the FormData
    data.append("publicationUrls", publicationsJsonData);

    const result: IStoreAndVerifyVusFileResponse = await fetch(
      `/vus/add_publications/${input.variantId}`,
      {
        method: "POST",
        body: data,
        // headers: {
        // "Content-Type": "multipart/form-data",
        // },
        cache: "no-store", //TODO: is it needed?
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async deleteVariant(
    input: IDeleteVariantRequest
  ): Promise<IDeleteVariantResponse> {
    const result: IDeleteVariantResponse = await fetch(
      `/vus/delete/${input.variantId}`,
      {
        method: "DELETE",
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

  async addSamples(input: IAddSamplesRequest): Promise<IAddSamplesResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("samplesToAdd", JSON.stringify(input.samplesToAdd));

    const result: IAddSamplesResponse = await fetch(
      `/vus/add-existing-samples/${input.variantId}`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async addNewSample(input: IAddNewSampleRequest): Promise<IAddNewSampleResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleToAdd", JSON.stringify(input.sampleToAdd));

    const result: IAddNewSampleResponse = await fetch(
      `/vus/add-new-sample/${input.variantId}`,
      {
        method: "POST",
        body: data,
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
