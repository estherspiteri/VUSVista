import customFetch from "../api/api.service";
import {
  IAddNewSampleRequest,
  IAddNewSampleResponse,
  IAddPublicationsRequest,
  IAddPublicationsResponse,
  IAddSamplesRequest,
  IAddSamplesResponse,
  ICheckFileForMultipleGenesRequest,
  ICheckFileForMultipleGenesResponse,
  ICheckFileForValidGenesRequest,
  ICheckFileForValidGenesResponse,
  ICheckFileUploadStatusesRequest,
  ICheckFileUploadStatusesResponse,
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
  IRemoveSamplesRequest,
  IRemoveSamplesResponse,
  IStoreAndVerifyVusFileRequest,
  IStoreAndVerifyVusFileResponse,
  IUpdateRsidRequest,
  IUpdateRsidResponse,
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

    const result: IUploadVusResponse = await customFetch(`/vus/upload`, {
      method: "POST",
      body: data,
      cache: "no-store",
    })
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error));

    return result;
  }

  async checkFileForMultipleGenes(
    input: ICheckFileForMultipleGenesRequest
  ): Promise<ICheckFileForMultipleGenesResponse> {
    let data = new FormData();
    data.append("file", input.vusFile);

    const result: ICheckFileForMultipleGenesResponse = await customFetch(
      `/vus/file/multiple-genes-check`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async checkFileForValidGenes(
    input: ICheckFileForValidGenesRequest
  ): Promise<ICheckFileForValidGenesResponse> {
    let data = new FormData();
    data.append("file", input.vusFile);

    const multipleGenesSelectionJsonData = input.multipleGenesSelection
      ? JSON.stringify(input.multipleGenesSelection)
      : "";

    // Append the JSON string as a blob to the FormData
    data.append("multipleGenesSelection", multipleGenesSelectionJsonData);

    const result: ICheckFileForMultipleGenesResponse = await customFetch(
      `/vus/file/existing-genes-check`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response) => {
        return response;
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

    const genesNotFoundSelectionJsonData = input.genesNotFoundSelection
      ? JSON.stringify(input.genesNotFoundSelection)
      : "";

    // Append the JSON string as a blob to the FormData
    data.append("genesNotFoundSelection", genesNotFoundSelectionJsonData);

    const result: IStoreAndVerifyVusFileResponse = await customFetch(
      `/vus/file`,
      {
        method: "POST",
        body: data,
        cache: "no-store",
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async checkFileUploadStatuses(
    input: ICheckFileUploadStatusesRequest
  ): Promise<ICheckFileUploadStatusesResponse> {
    const result: ICheckFileUploadStatusesResponse = await customFetch(
      `/vus/file/check-status/${input.taskIds.join(",")}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
        cache: "no-store",
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async loadAllVus(): Promise<ILoadAllVusResponse> {
    const result: ILoadAllVusResponse = await customFetch(`/vus/view`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
      cache: "no-store",
    })
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getVus(input: IGetVusRequest): Promise<IGetVusResponse> {
    const result: IGetVusResponse = await customFetch(
      `/vus/view/${input.vusId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
        cache: "no-cache",
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async verifyGene(input: IVerifyGeneRequest): Promise<IVerifyGeneResponse> {
    const result: IVerifyGeneResponse = await customFetch(
      `/vus/gene/${input.geneName}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getAllAcmgRules(): Promise<IGetAllAcmgRulesResponse> {
    const result: IGetAllAcmgRulesResponse = await customFetch(
      `/vus/all-acmg-rules`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getClinvarUpdates(
    input: IGetClinvarUpdatesRequest
  ): Promise<IGetClinvarUpdatesResponse> {
    const result: IGetClinvarUpdatesResponse = await customFetch(
      `/vus/get_clinvar_updates/${input.clinvarId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
        cache: "no-cache",
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getPublicationUpdates(
    input: IGetPublicationUpdatesRequest
  ): Promise<IGetPublicationUpdatesResponse> {
    const result: IGetPublicationUpdatesResponse = await customFetch(
      `/vus/get_publication_updates/${input.variantId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
        cache: "no-cache",
      }
    )
      .then((response) => {
        return response;
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

    const result: IStoreAndVerifyVusFileResponse = await customFetch(
      `/vus/add_publications/${input.variantId}`,
      {
        method: "POST",
        body: data,
        cache: "no-store",
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async deleteVariant(
    input: IDeleteVariantRequest
  ): Promise<IDeleteVariantResponse> {
    const result: IDeleteVariantResponse = await customFetch(
      `/vus/delete/${input.variantId}`,
      {
        method: "DELETE",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async addSamples(input: IAddSamplesRequest): Promise<IAddSamplesResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("samplesToAdd", JSON.stringify(input.samplesToAdd));

    const result: IAddSamplesResponse = await customFetch(
      `/vus/add-existing-samples/${input.variantId}`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async addNewSample(
    input: IAddNewSampleRequest
  ): Promise<IAddNewSampleResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleToAdd", JSON.stringify(input.sampleToAdd));

    const result: IAddNewSampleResponse = await customFetch(
      `/vus/add-new-sample/${input.variantId}`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async removeSamples(
    input: IRemoveSamplesRequest
  ): Promise<IRemoveSamplesResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleIdsToRemove", JSON.stringify(input.sampleIdsToRemove));

    const result: IRemoveSamplesResponse = await customFetch(
      `/vus/remove-samples/${input.variantId}`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async updateRsid(input: IUpdateRsidRequest): Promise<IUpdateRsidResponse> {
    const result: IUpdateRsidResponse = await customFetch(
      `/vus/update-rsid/${input.variantId}/${input.newRsid}`,
      {
        method: "POST",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const vusService = new VusService();
